using ModelingToolkit, Plots, DifferentialEquations, LinearAlgebra
const Φ₀ = 2.067833848e-15 #flux quantum


@variables t #time variable
const D = Differential(t) #define an operator for the differentiation w.r.t. time

function build_component(;name) #general component with phase and current states
    sts = @variables θ(t)=1.0 i(t)=1.0 iₘ(t) = 1.0          #iₘ is global, not unique to components
    ODESystem(Equation[], t, sts, []; name=name)
end

struct Component
    sys::ODESystem
end

# σB given as an array [0 0 -1] representing which loop delivers current + direction
function build_resistor(;name, R = 1.0, σB = []) #Builds ODESystem for resistor using Component
    @named component = build_component()
    @unpack θ, i, iₘ = component
    ps = @parameters R=R σB=σB
    iₘs =  0
    for i in 1:length(σB)
        iₘs = iₘs + σB[i]*iₘ[i]
    end
    eqs = [
            #i~*D(θ)*Φ₀/(2*pi*R)
            D(θ) ~ i*(2*pi*R)/Φ₀
            i ~ iₘs
          ]
    sys = extend(ODESystem(eqs, t, [], ps; name=name), component)
    Component(sys)
end

function build_capacitor(;name, C = 1.0, σB = []) #builds ODESystem for capacitor using Component
    @named component = build_component()
    @unpack θ, i, iₘ = component
    ps = @parameters C=C σB=σB
    iₘs =  0
    for i in 1:length(σB)
        iₘs = iₘs + σB[i]*iₘ[i]
    end
    eqs = [
            #i~D(D(θ))*(Φ₀*C)/(2*pi)
            D(D(θ))~i*2*pi/(Φ₀*C)
            i ~ iₘs
          ]
    sys = extend(ODESystem(eqs, t, [], ps; name=name), component)
    Component(sys)
end

function build_JJ(;name, Io = 1.0, R = 1.0, C = 1.0, σB = []) #Builds ODESystem for resistor using Component
    @named component = build_component()
    @unpack θ, i, iₘ = component
    ps = @parameters Io=Io R=R C=C σB=σB
    iₘs =  0
    for i in 1:length(σB)
        iₘs = iₘs + σB[i]*iₘ[i]
    end
    eqs = [
            i ~ Io*sin(θ) + D(θ)*Φ₀/(2*pi*R) + D(D(θ))*(Φ₀*C)/(2*pi)
            #D(θ)~i*(2*pi*R)/Φ₀
            i ~ iₘs
          ]
    sys = extend(ODESystem(eqs, t, [], ps; name=name), component)
    Component(sys)
end

struct Loop
    sys::ODESystem
end

function build_loop(;name, Φₑ = 0.0, Lᵢ = [] σA = [])   # NUMBER OF COMPONENTS NEEDED
    sts = @variables iₘ(t) = 1.0 θ(t)=1.0               # θ[1:3](t)
    ps = @parameters Φₑ=Φₑ Lᵢ=Lᵢ σA=σA
    θs = 0
    iₘs = 0
    for i in 1:length(σA)
        θs = θs + σA[i]*θ[i]
        iₘs = iₘs + Lᵢ[i]*iₘ[i]
    end
    eqs = [
            θs ~ Φₑ - iₘs
          ]
    sys = ODESystem(eqs, t, sts, ps, name=name)
    Loop(sys)
end

compose(everthing)




#=
struct CurrentSource
    sys::ODESystem
end

function build_current_loop(;name, I = 1.0)
    @named loop = build_loop()
    @unpack iₘ, Φₗ = loop
    ps = @parameters I = I
    eqs = [I~iₘ]                # easier to define these for each component???  ~~~
    sys = extend(ODESystem(eqs,t, [], ps, name=name), loop)
    CurrentSource(sys)
end

declared_components = []
branched_components = []

function loop_equations(l::Loop, σ::Vector{Int}, cs::Component ...)
    cnames = [c for c in map(c->nameof(c.sys),cs)]
       for c in cnames
           if c in declared_components
               continue
           else
               push!(declared_components,c)
           end
       end
    θ_eqs = [dot([c.sys.θ for c in cs],σ) ~ l.sys.Φₑ-l.sys.Φₗ] #im needs to be multiplied? ~~~
    i_eqs = map((c,Σ)->Σ*c.sys.i~l.sys.iₘ,cs,σ)

    return append!(θ_eqs, i_eqs)
end




function connect_loops(l1::Loop, l2::Loop, σ::Vector{Int}, cs::Component ...; Lᵢⱼ = 0.0)
    cnames = [c for c in map(c->nameof(c.sys),cs)]
       for c in cnames
           if c in branched_components
               continue
           else
               push!(branched_components,c)
           end
       end
    ps = @parameters Lᵢⱼ = Lᵢⱼ
    loop_eqs = [
    l1.sys.Φₗ ~ -l1.sys.Lᵢᵢ * l1.sys.iₘ + Lᵢⱼ * l2.sys.iₘ,      #These are already known in L matrix ~~~
    l2.sys.Φₗ ~ -l2.sys.Lᵢᵢ * l2.sys.iₘ + Lᵢⱼ * l1.sys.iₘ
    ]



    return append!(loop_eqs,branch_eqs), ps
end
=#