using JLD2, FileIO, ModelingToolkit, Plots, DifferentialEquations, LinearAlgebra, Statistics
include("model_builder.jl")

function open_file(name)
    global f_name = name
    try
        if endswith(f_name, ".jld2")
            file = jldopen(f_name, "r")
        else 
            file = jldopen("$f_name.jld2", "r")
        end
    catch e
        if isa(e, SystemError)
            println("File does not exist")
        end
        #break?
    end
    if endswith(f_name, ".jld2")
        file = jldopen(f_name, "r")
    else 
        file = jldopen("$f_name.jld2", "r")
    end
    #file = jldopen("$f_name.jld2", "r")
    global numLoops = read(file, "editing/numLoops")
    global componentLoopDict = read(file, "editing/componentLoopDict")
    global CPD = read(file, "editing/componentParamDict")
    global mutualInd = read(file, "editing/mutualInd")
    global junctions = read(file, "editing/junctions")
    global loops = read(file, "editing/loops")
    global k = read(file, "matrices/k")
    global L = read(file, "matrices/L")
    global σA = read(file, "matrices/σA")
    global σB = read(file, "matrices/σB")
    global componentPhaseDirection = read(file, "matrices/componentPhaseDirection")
    close(file)
end

### Ask user to input any variables that are still in symbolic form
function symbolic_assign()
    symbolDict = Dict()
    for comp in CPD
        if (comp[1][1] == 'J')
            for i in 1:3
                if isa(comp[2][i], Symbol)
                    symbolDict[comp[2][i]] = push!(get(symbolDict, comp[2][i], []), (comp[1], i))
                end
            end
        else
            if isa(comp[2], Symbol)
                symbolDict[comp[2]] = push!(get(symbolDict, comp[2], []), comp[1])
            end
        end
    end
    for sym in symbolDict
        println("Please enter a value for $(sym[1])")
        input = readline()
        for i in sym[2]
            if isa(i, Tuple)
                CPD[i[1]][i[2]] = parse(Float64, input)
            else
                CPD[i] = input
            end
        end
    end
end

#Build the circuit based on the open file
function build()
    eqs = Equation[]

    built_loops = []
    external_flux_strength = Φ₀/10        #### THIS NEEDS TO BE SET BY USER

    for i in 1:numLoops
        println("loop $(i-1)")
        if (startswith(loops[i][1], "I"))
            new_l = "@named loop$(i-1) = build_current_source_loop(I = $(CPD["Ib"]))"
            new_l = Meta.parse(new_l)
            new_l = eval(new_l)
        else
            new_l = "@named loop$(i-1) = build_loop(Φₑ = $(external_flux_strength*k[i]))"
            new_l = Meta.parse(new_l)
            new_l = eval(new_l)
        end
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

    old_sys = []
    global u0 = Pair{Num, Float64}[]

    for comp in built_components
        push!(old_sys, comp[2].sys)
        if (comp[1][1] != 'V')
            push!(u0, comp[2].sys.θ=>0.0)
            push!(u0, comp[2].sys.i=>0.0)
            if (uppercase(comp[1][1]) == 'C')
                push!(u0, D(comp[2].sys.θ)=>0.0)
            elseif (uppercase(comp[1][1]) == 'J')
                push!(u0, D(comp[2].sys.θ)=>0.0)
            end
        end
    end
    for loop in built_loops
        push!(old_sys, loop.sys)
    end
    sys = Vector{ODESystem}(old_sys)

    inductance(eqs, L, built_loops)
    current_flow(eqs, componentPhaseDirection, built_loops, built_components)
    for i in 1:numLoops
        add_loop!(eqs, built_loops[i], σA[i,:], built_components)
    end

    @named _model  =  ODESystem(eqs, t)

    @named model = compose(_model, sys)

    println()
    display(equations(model))
    println()
    display(states(model))
    println()
    global new_model = structural_simplify(model)
end

#Solve initial conditions
function solve_init()
    tspan_ini = (0.0,1e-9)
    prob = ODEProblem(new_model, u0, tspan_ini, save_everystep = false, progress=true)
    sol = solve(prob, ROS3P())
    global u0 = sol[:,end]
end

#Give a timespan for the simulation and solve
function tspan(t_init, t_end)
    tspan = (parse(Float64, t_init), parse(Float64, t_end))
    tsaves = LinRange(tspan[1],tspan[2], 50000)
    global prob = ODEProblem(new_model, u0, tspan, saveat=tsaves, progress=true)
    global sol = solve(prob)
end

#Plot a single component
function single_plot(comp, param)
    try
        if (param == "i")
            str = "$comp.sys.i"
            ex = Meta.parse(str)
            p = plot(sol, vars=[eval(ex)])
            png("$(f_name)_$(comp)_i")
        elseif (param == "V")
            str = "D($comp.sys.θ)"
            ex = Meta.parse(str)
            p = plot(sol, vars=[eval(ex)])
            png("$(f_name)_$(comp)_V")
        end
    catch e
        if isa(e, UndefVarError)
            println(" --- Component does not exist ---")
        elseif isa(e, ArgumentError)
            println(" --- Cannot plot voltage through resistors at this point in time ---")
        end
    end
end

#Ensemble functions, not complete
function prob_func(prob,i,repeat) #problem funtion modifies input parameters
    global flux_vec = LinRange(0.0,10,10)
    new_p = prob.p[1:numLoops]*flux_vec[i]
    remake(prob,p=[new_p; prob.p[(numLoops+1):end]])
end
function output_func(sol, i)
    (mean(Φ₀/(2.0*pi)*1.0e+6*sol),false)
end# output_func

#=function ensemble()
    ensemble_prob = EnsembleProblem(prob, prob_func=prob_func); #ensemble

    sim = solve(ensemble_prob,Tsit5(),EnsembleDistributed(),trajectories=200)
    #sim = solve(ensemble_prob,Vern6(),EnsembleThreads(),trajectories=N_fvalues, maxiters=1e9, progress=true) #solve ensemble
    plot(sim)
end=#

#=  Commands running from Julia REPL
open_file("apf")
build()
solve_init()
tspan("0.0", "1e-12")
single_plot("Ra","i")

ensemble_prob = EnsembleProblem(prob, prob_func=prob_func)#, output_func = output_func); #ensemble
sim = solve(ensemble_prob,EnsembleThreads(),trajectories=10, maxiters=1e9, progress=true) #solve ensemble
plot(sim, vars=(D(J1.sys.θ), J1.sys.i))


graph = plot(flux_vec,sim[:],
    title = "RLC-SQUID Parallel Config \n Flux-Voltage Response",
    xlabel = "Flux (Φe/Φ₀)",
    ylabel = "Averaged Vout (μVrms)"
    )
display(graph)

plot(sim)
#ensemble()

#= Run from terminal commands
while true
    println(" --- Enter 'help' if you need help --- ")
    input = readline()
    if (input == "~")
        break
    elseif (lowercase(input) == "help")
        #Give instructions on how to use
    elseif startswith(lowercase(input), "open_file")
        open_file(input[11:end-1])  
    elseif startswith(lowercase(input), "build")
        build()
    elseif startswith(lowercase(input), "solve_init")
        solve_init()
    elseif startswith(lowercase(input), "tspan")
        ts = split(strip(input[7:end-1]),',')   
        tspan(strip(ts[1]), strip(ts[2]))               #strip() removes spaces
    elseif startswith(lowercase(input), "single_plot")
        s_plot = split(input[13:end-1], ',')
        single_plot(strip(s_plot[1]), strip(s_plot[2])) #strip() removes spaces
    end
end
=#
 
prob.p
prob.p[1:numLoops]
new_p = prob.p[1:numLoops]*2
newer_p = [new_p; prob.p[4:end]]
