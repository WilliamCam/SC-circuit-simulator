using ModelingToolkit, Plots, DifferentialEquations, LinearAlgebra
const Φ₀ = 2.067833848e-15 #flux quantum

@variables t #time variable
D = Differential(t)
D2 = Differential(t)^2

function build_component(;name) #general component with phase and current states
    sts = @variables θ(t)=1.0 i(t)=1.0
    ODESystem(Equation[], t, sts, []; name=name)
end

struct Component
    sys::ODESystem
end

function build_resistor(;name, R = 1.0) #Builds ODESystem for resistor using Component
    @named component = build_component()
    @unpack θ, i = component
    ps = @parameters R=R
    eqs = [
            D(θ)~i*(2*pi*R)/Φ₀
          ]
    sys = extend(ODESystem(eqs, t, [], ps; name=name), component)
    Component(sys)
end

#=function build_port(;name, R = 10.0) #Builds ODESystem for resistor using Component
    @named component = build_component()
    @unpack θ, i = component
    ps = @parameters R=R
    eqs = [
            D(θ)~i*(2*pi*R)/Φ₀
          ]
    sys = extend(ODESystem(eqs, t, [], ps; name=name), component)
    Component(sys)
end=#

function build_capacitor(;name, C = 1.0) #builds ODESystem for capacitor using Component
    @named component = build_component()
    @unpack θ, i = component
    ps = @parameters C=C
    eqs = [
            D2(θ)~i*2*pi/(Φ₀*C)
          ]
    sys = extend(ODESystem(eqs, t, [], ps; name=name), component)
    Component(ode_order_lowering(sys))
end

function build_JJ(;name, I0 = 1.0, R = 1.0, C = 1.0) #builds ODESystem for JosephsonJunction using Component
    @named component = build_component()
    @unpack θ, i = component
    ps = @parameters C=C
    eqs = [
            D2(θ) ~ (i - I0*sin(θ) - D(θ)*Φ₀/(2*pi*R))*(2*pi)/(Φ₀*C)
          ]
    sys = extend(ODESystem(eqs, t, [], ps; name=name), component)
    Component(ode_order_lowering(sys))
end

function build_voltage_source(;name, V = 1.0, ω = 0.0)
    @named component = build_component()
    @unpack θ, i = component
    ps = @parameters V=V
    eqs = [
            D(θ)~ - V*cos(ω*t)*2*pi/Φ₀
          ]
    sys = extend(ODESystem(eqs, t, [], ps; name=name), component)
    Component(sys)
end

#=function build_drive_port(;name, V = 1.0, ω = 0.0)
    @named component = build_component()
    @unpack θ, i = component
    ps = @parameters V=V
    eqs = [
            D(θ)~ - V*cos(ω*t)*2*pi/Φ₀
          ]
    sys = extend(ODESystem(eqs, t, [], ps; name=name), component)
    Component(sys)
end=# #How to include eqs for resistor?

struct Loop
    sys::ODESystem
end

function build_loop(;name, Φₑ = 0.0)
    sts = @variables iₘ(t) = 0.0 Φₗ(t)=0.0
    ps = @parameters Φₑ = Φₑ
    sys = ODESystem(Equation[], t, sts, ps, name=name)
    Loop(sys)
end

struct CurrentSourceLoop
    sys::ODESystem
end

function build_current_source_loop(;name, I = 1.0, ω = 0.0)
    @named loop = build_loop()
    @unpack sys = loop
    @unpack iₘ, Φₗ = sys
    ps = @parameters I = I
    eqs = [
            0 ~ iₘ - I*cos(ω*t)
          ]
    sys = extend(ODESystem(eqs, t, [], ps, name=name), loop.sys)
    CurrentSourceLoop(sys)
end

struct ComponentFlow
    c::Component
    σ::Float64
end

function add_loop!(
    eqs::Vector{Equation},l, σ::Vector, cs
    )
    push!(eqs,0 ~ l.sys.Φₑ-l.sys.Φₗ - dot([Φ₀/(2*pi)*c.sys.θ for (n, c) in cs], σ)) #<---- σA . θ
end

function inductance(eqs::Vector{Equation}, L, built_loops)
    for i in 1:length(built_loops)                              #Starts at loop1 as loop0 (Ib) can not have inductance
        if (string(built_loops[i])[1:4] == "Loop")
            push!(eqs, built_loops[i].sys.Φₗ ~ dot(L[:,i], [l.sys.iₘ for l in built_loops]))
        end
    end
end

function current_flow(eqs::Vector{Equation}, componentPhaseDirection, built_loops, built_components)
    for (comp, phase_array) in componentPhaseDirection              #Loop through each junction where "comp" is the junction name and "phase_array" is the row of σB corresponding to the currente junction
        σB_im = 0                                                   #Temp variable to store σB * im for the current junction
        for i in 1:length(phase_array)
            σB_im = σB_im + phase_array[i]*built_loops[i].sys.iₘ    #Adds direction of im from each loop through the junction
        end
        push!(eqs, 0 ~ σB_im - get(built_components, comp, -1).sys.i)   #Current of junction i = σB_im (direction of current from each loop through junction)
    end
end