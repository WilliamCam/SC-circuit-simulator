using JLD2, FileIO, ModelingToolkit, Plots, DifferentialEquations, LinearAlgebra
include("model_builder.jl")

f_name = "RC.jld2"
file = jldopen(f_name, "r")
loops = load(f_name, "editing/loops")
CPD = load(f_name, "editing/componentParamDict")
CLD = load(f_name, "editing/componentLoopDict")
junctions = load(f_name, "editing/junctions")
numLoops = load(f_name, "editing/numLoops")
mutualInd = load(f_name, "editing/mutualInd")
k = load(f_name, "matrices/k")
σB = load(f_name, "matrices/σB")
componentPhaseDirection = load(f_name, "matrices/componentPhaseDirection")
σA = load(f_name, "matrices/σA")
L = load(f_name, "matrices/L")
close(file)

j_len = length(junctions)

tspan = (0.0,100.0)
eqs = Equation[]

built_loops = []    #build loops
for i in 1:numLoops
    @named loop = build_loop()
    push!(built_loops, loop)
end

built_components = Dict()        #Store built components
for j in junctions
    println(j)
    if (j[1] == 'R')
        param = get(CPD, j, 0)
        new_c = "@named $j = build_resistor(R = $param)"
        new_c = Meta.parse(new_c)
        new_c = eval(new_c)
        built_components[j] = new_c
    elseif (j[1] == 'C')
        param = get(CPD, j, 0)
        new_c = "@named $j = build_capacitor(C = $param)"
        new_c = Meta.parse(new_c)
        new_c = eval(new_c)
        built_components[j] = new_c
    elseif (j[1] == 'V')
        param = get(CPD, j, 0)
        new_c = "@named $j = build_voltage_source(V = $param)"
        new_c = Meta.parse(new_c)
        new_c = eval(new_c)
        built_components[j] = new_c
    end
end

loops_built_components = []
for i in 1:numLoops
    current_components = Component[]
    for j in loops[i]
        push!(current_components, get(built_components, j, 0))
    end
    push!(loops_built_components, current_components)
end

D = Differential(t)

old_sys = []
my_u0 = Pair{Num, Float64}[]
for comp in built_components
    push!(old_sys, comp[2].sys)
    if (comp[1][1] != 'V')
        push!(my_u0, comp[2].sys.θ=>0.0)
        push!(my_u0, comp[2].sys.i=>0.0)
        if (comp[1][1] == 'C')
            push!(my_u0, D(comp[2].sys.θ)=>0.0)
        end
    end
end
for loop in built_loops
    push!(old_sys, loop.sys)
end

sys = Vector{ODESystem}(old_sys)

current_flow(eqs, componentPhaseDirection, built_loops, built_components)
for i in 1:numLoops
    add_loop!(eqs, built_loops[i], σA[i,:], loops_built_components[i])
end

@named _model  =  ODESystem(eqs, t)

@named model = compose(_model, sys)

new_model = structural_simplify(model)

new_new_model = ode_order_lowering(new_model)

tsaves = LinRange(tspan[1],tspan[2], Integer(tspan[2]*50))

prob = ODEProblem(new_new_model, my_u0, tspan, saveat=tsaves)

sol = solve(prob, Rodas4(), maxiters=1e6, abstol = 1e-9, reltol=1e-12)

plot(sol, vars=[R1.sys.i])
plot(sol.t, (Φ₀/2*pi)*sol[D(C1.sys.θ)])
