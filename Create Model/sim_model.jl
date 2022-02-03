using JLD2, FileIO, ModelingToolkit, Plots, DifferentialEquations, LinearAlgebra
include("model_builder.jl")

f_name = "vsource_jj.jld2"
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

CPD["J1"] = [1.0, 10, 3.2910597840193506e-14, 1]
CPD["V1"] = 1.0e-6
#phi0 = 2.067833848e-15
#I₀₁ = 10.0e-6
#Gj1 = 1/10
#βc₁ = 0.1
#Cj1 = βc₁*phi0/(2*pi)*Gj1^2/(I₀₁)

j_len = length(junctions)

eqs = Equation[]

built_loops = []    #build loops
for i in 1:numLoops
    new_l = "@named loop$(i-1) = build_loop()"
    new_l = Meta.parse(new_l)
    new_l = eval(new_l)
    push!(built_loops, new_l)
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
    elseif (j[1] == 'J')
        params = get(CPD, j, 0)
        new_c = "@named $j = build_JJ(I0 = $(params[1]), R = $(params[2]), C = $(params[3]))"
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

#loops_built_components[2]

D = Differential(t)

old_sys = []
my_u0 = Pair{Num, Float64}[]
for comp in built_components
    push!(old_sys, comp[2].sys)
    if (comp[1][1] != 'V')
        push!(my_u0, comp[2].sys.θ=>0.0)
        push!(my_u0, comp[2].sys.i=>0.0)
        if (uppercase(comp[1][1]) == 'C')
            push!(my_u0, D(comp[2].sys.θ)=>0.0)
        elseif (uppercase(comp[1][1]) == 'J')
            push!(my_u0, D(comp[2].sys.θ)=>0.0)
        end
    end
end
for loop in built_loops
    push!(old_sys, loop.sys)
end

my_u0

sys = Vector{ODESystem}(old_sys)

#σA
#σA[1,:]
#σA[2,:]

#built_components
#for (n, c) in built_components
#    println(n)
#end

current_flow(eqs, componentPhaseDirection, built_loops, built_components)
for i in 1:numLoops
    add_loop!(eqs, built_loops[i], σA[i,:], built_components)
end

@named _model  =  ODESystem(eqs, t)

@named model = compose(_model, sys)

new_model = structural_simplify(model)


for eq in new_new_model.eqs
    println()
    println(eq)
end

tspan = (0.0,1e-9)
tsaves = LinRange(tspan[1],tspan[2], 5000)
prob = ODEProblem(new_model, my_u0, tspan, saveat=tsaves)
sol = solve(prob, ROS3P())#, maxiters=1e6, abstol = 1e-9, reltol=1e-12)

plot(sol, vars=[J1.sys.i])
#plot(sol, vars=[D(J1.sys.θ)])
#plot(sol.t, (Φ₀/2*pi)*sol[D(C1.sys.θ)])
