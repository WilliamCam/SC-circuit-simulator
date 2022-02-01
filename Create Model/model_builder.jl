using ModelingToolkit, Plots, DifferentialEquations, LinearAlgebra
const Φ₀ = 2.067833848e-15 #flux quantum


@variables t #time variable
D = Differential(t)
D2 = Differential(t)^2

function build_component(;name) #general component with phase and current states
    sts = @variables θ(t)=1.0 i(t)=1.0
    ODESystem(Equation[], t, sts, []; name=name)
end

function build_JJ_component(;name) #general component with phase and current states
    sts = @variables θ(t)=1.0 i(t)=1.0 RHS(t)=1.0
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
            #i~*D(θ)*Φ₀/(2*pi*R)
            D(θ)~i*(2*pi*R)/Φ₀,
          ]
    sys = extend(ODESystem(eqs, t, [], ps; name=name), component)
    Component(sys)
end

function build_capacitor(;name, C = 1.0) #builds ODESystem for capacitor using Component
    @named component = build_component()
    @unpack θ, i = component
    ps = @parameters C=C
    eqs = [
            #i~D(D(θ))*(Φ₀*C)/(2*pi)
            D2(θ)~i*2*pi/(Φ₀*C)
          ]
    sys = extend(ODESystem(eqs, t, [], ps; name=name), component)
    Component(sys)
end

function build_JJ(;name, Io = 1.0, R = 1.0, C = 1.0) #builds ODESystem for JosephsonJunction using Component
    @named component = build_JJ_component()
    @unpack θ, i, RHS = component
    ps = @parameters C=C
    eqs = [
            #i~D(D(θ))*(Φ₀*C)/(2*pi)
            D(θ) ~ (i - Io*sin(θ) - RHS)*(2*pi*R)/Φ₀
            D2(θ) ~ RHS*(2*pi)/(Φ₀*C)
          ]
    sys = extend(ODESystem(eqs, t, [], ps; name=name), component)
    Component(sys)
end

function build_voltage_source(;name, V = 1.0)
    @named component = build_component()
    @unpack θ, i = component
    ps = @parameters V=V
    eqs = [
            D(θ)~ - V*2*pi/Φ₀
          ]
    sys = extend(ODESystem(eqs, t, [], ps; name=name), component)
    Component(sys)
end

struct Loop
    sys::ODESystem
end

function build_loop(;name, Φₑ = 0.0, Lᵢᵢ = 0.0)
    sts = @variables iₘ(t) = 0.0 Φₗ(t)=0.0
    ps = @parameters Φₑ = Φₑ Lᵢᵢ = Lᵢᵢ
    sys = ODESystem(Equation[], t, sts, ps, name=name)
    Loop(sys)
end

struct CurrentSourceLoop
    sys::ODESystem
end

function build_current_source_loop(;name, I = 1.0, Lᵢᵢ = 0.0)
    @named loop = build_loop()
    @unpack sys = loop
    @unpack iₘ, Φₗ, Lᵢᵢ = sys
    ps = @parameters I = I
    eqs = [I~iₘ]
    sys = extend(ODESystem(eqs,t, [], ps, name=name), loop.sys)
    CurrentSourceLoop(sys)
end

declared_cs = []

struct ComponentFlow
    c::Component
    σ::Float64
end

function add_loop!(
    eqs::Vector{Equation},l::Loop, σ::Vector, cs
    )
    #=branched_cs = []
    names = [cname for cname in map(c->nameof(c.sys),cs)]
    new_cs = map((c,σⱼ)->ComponentFlow(c,σⱼ),cs,σ)
       for cFlow in new_cs
           if nameof(cFlow.c.sys) in declared_cs
               push!(branched_cs,cFlow)
           else
               push!(declared_cs,nameof(cFlow.c.sys))
           end
       end=##

    push!(eqs,dot([Φ₀/(2*pi)*c.sys.θ for (n, c) in cs],σ) ~ l.sys.Φₑ-l.sys.Φₗ)

    #=for cFlow in setdiff(new_cs,branched_cs)
        push!(eqs, cFlow.σ*cFlow.c.sys.i~l.sys.iₘ)
    end

    for cFlow in branched_cs
        eqs = substitute(eqs,Dict(cFlow.c.sys.i=>cFlow.c.sys.i-cFlow.σ*l.sys.iₘ))
    end=#
    #^^^ ---- Form these eqs from componentPhaseDirection 

    push!(eqs, l.sys.Φₗ ~ -l.sys.Lᵢᵢ*l.sys.iₘ)
end

function add_current_source_loop!(
    eqs::Vector{Equation}, lc::CurrentSourceLoop, σ::Vector{Float64}, cs::Component ...
    )
    branched_cs = []
    names = [cname for cname in map(c->nameof(c.sys),cs)]
    new_cs = map((c,σⱼ)->ComponentFlow(c,σⱼ),cs,σ)
       for cFlow in new_cs
           if nameof(cFlow.c.sys) in declared_cs
               push!(branched_cs,cFlow)
           else
               push!(declared_cs,nameof(cFlow.c.sys))
           end
       end


       for cFlow in setdiff(new_cs,branched_cs)
           push!(eqs, cFlow.σ*cFlow.c.sys.i~lc.sys.iₘ)
       end

       for cFlow in branched_cs
           eqs = substitute(eqs,Dict(cFlow.c.sys.i=>cFlow.c.sys.i-cFlow.σ*lc.sys.iₘ))
       end

       push!(eqs, lc.sys.Φₗ ~ -lc.sys.Lᵢᵢ*lc.sys.iₘ)
end

function mutual_inductance(eqs, l1::Loop, l2::Loop; Lᵢⱼ::Float64 = 1.0)
    substitute(eqs, Dict(-l1.sys.Lᵢᵢ*l1.sys.iₘ=>-l1.sys.Lᵢᵢ*l1.sys.iₘ-l2.sys.iₘ*Lᵢⱼ))
    substitute(eqs, Dict(-l2.sys.Lᵢᵢ*l2.sys.iₘ=>-l2.sys.Lᵢᵢ*l2.sys.iₘ-l1.sys.iₘ*Lᵢⱼ))
end

function current_flow(eqs::Vector{Equation}, componentPhaseDirection, built_loops, built_components)
    for (comp, phase_array) in componentPhaseDirection              #Loop through each junction where "comp" is the junction name and "phase_array" is the row of σB corresponding to the currente junction 
        σB_im = 0                                                   #Temp variable to store σB * im for the current junction
        for i in 1:length(phase_array)
            σB_im = σB_im + phase_array[i]*built_loops[i].sys.iₘ    #Adds direction of im from each loop through the junction
        end
        push!(eqs, get(built_components, comp, -1).sys.i ~ σB_im)   #Current of junction i = σB_im (direction of current from each loop through junction)
    end
end