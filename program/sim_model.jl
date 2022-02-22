using JLD2, FileIO, ModelingToolkit, Plots, DifferentialEquations, LinearAlgebra, Statistics, Dates
include("model_builder.jl")

function help()
    println("This is helpful")
end

function open_file(name)
    f_name = name
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
    numLoops = read(file, "editing/numLoops")
    #componentLoopDict = read(file, "editing/componentLoopDict")
    CPD = read(file, "editing/componentParamDict")
    #mutualInd = read(file, "editing/mutualInd")
    junctions = read(file, "editing/junctions")
    loops = read(file, "editing/loops")
    k = read(file, "matrices/k")
    L = read(file, "matrices/L")
    σA = read(file, "matrices/σA")
    #σB = read(file, "matrices/σB")
    componentPhaseDirection = read(file, "matrices/componentPhaseDirection")
    close(file)
    return numLoops, CPD, junctions, loops, k, L, σA, componentPhaseDirection
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
        println("Symbol: \n$(sym[1])\nComponents")
        display(sym[2])
        println("Please enter a value for $(sym[1])")
        input = readline()
        for i in sym[2]
            if isa(i, Tuple)
                CPD[i[1]][i[2]] = parse(Float64, input)
            else
                CPD[i] = parse(Float64, input)
            end
        end
    end
end

#Edit parameters while in program - Does not work when using Julia Repl
function edit_param()
    while true
        display(CPD)
        println()
        println("Enter a component followed by '-' followed by it's new parameter. E.g. to change Ra from 2 Ohm to 3 Ohm enter Ra-3\nEnter ~ when finished editing parameters")
        input = readline()
        if (input == "~")
            break
        end
        comp_param = split(input, '-')
        if (comp_param[1] in keys(CPD))
            CPD[comp_param[1]] = comp_param[2]
        else
            println("Component does not exist, try agian")
        end
    end
end

#Build the circuit based on the open file
function build(Φext)
    eqs = Equation[]

    built_loops = []

    for i in 1:numLoops
        println("loop $(i-1)")
        current_loop = false
        c_source = ""
        for comp in loops[i]
            if startswith(comp, 'I')
                current_loop = true
                c_source = comp
            end
        end
        if (current_loop)
            param = CPD[c_source]
            if (length(param) == 1)
                new_l = "@named loop$(i-1) = build_current_source_loop(I = $(param))"
            else 
                new_l = "@named loop$(i-1) = build_current_source_loop(I = $(param[1]), ω = $(1/param[2]))"
            end
            new_l = Meta.parse(new_l)
            new_l = eval(new_l)
        else
            new_l = "@named loop$(i-1) = build_loop(Φₑ = $(Φext*k[i]))"
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
            if (length(param) == 1)
                new_c = "@named $j = build_voltage_source(V = $param)"
            else
                new_c = "@named $j = build_voltage_source(V = $(param[1]), ω = $(1/param[2]))"
            end
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
    u0 = Pair{Num, Float64}[]

    for comp in built_components
        push!(old_sys, comp[2].sys)
        if (comp[1][1] != 'V')
            push!(u0, comp[2].sys.θ=>0.0)
            push!(u0, comp[2].sys.i=>0.0)
            if (uppercase(comp[1][1]) in ['C', 'J', 'V'])
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
    new_model = structural_simplify(model)
    return new_model, u0
end

#Solve initial conditions
function solve_init(old_u0)
    tspan_ini = (0.0,1e-9)
    prob = ODEProblem(new_model, old_u0, tspan_ini, save_everystep = false, progress=true)
    sol = solve(prob, ROS3P())
    new_u0 = sol[:,end]
    return new_u0
end

#Give a timespan for the simulation and solve
function tspan(t_init, t_end)
    tspan = (parse(Float64, t_init), parse(Float64, t_end))
    tsaves = LinRange(tspan[1],tspan[2], 50000)
    prob = ODEProblem(new_model, u0, tspan, saveat=tsaves, progress=true)
    sol = solve(prob)
    return sol
end

#Plot a single component
function single_plot(comp, param)
    try
        if (param == "i")
            str = "$comp.sys.i"
            ex = Meta.parse(str)
            p = plot(sol, vars=[eval(ex)])
            display(p)
            png("$(now())")
        elseif (param == "V")
            str = "D($comp.sys.θ)"
            ex = Meta.parse(str)
            p = plot(sol, vars=[eval(ex)])
            display(p)
            png("$(now())")
        end
    catch e
        if isa(e, UndefVarError)
            println(" --- Component does not exist ---")
        elseif isa(e, ArgumentError)
            dtime = sol.t[2]-sol.t[1]
            str = "$comp.sys.θ"
            ex = Meta.parse(str)
            dv = sol[eval(ex)]
            res = []
            for i in 1:length(dv)-1
                d = Φ₀/2pi*(dv[i+1] - dv[i])/dtime
                push!(res, d)
            end
            p = plot(sol.t[1:end-1], res)
            display(p)
            png("$(now())")
        end
    end
end


#=Ensemble functions, not complete
function prob_func(prob,i,repeat) #problem funtion modifies input parameters
    flux_vec = LinRange(0.0,10,10)
    new_p = prob.p[1:numLoops]*flux_vec[i]
    remake(prob,p=[new_p; prob.p[(numLoops+1):end]])
end
function output_func(sol, i)
    (mean(Φ₀/(2.0*pi)*1.0e+6*sol),false)
end# output_func
function ensemble()
    ensemble_prob = EnsembleProblem(prob, prob_func=prob_func, output_func = output_func) #ensemble
    sim = solve(ensemble_prob,Tsit5(),EnsembleDistributed(),trajectories=100)
    plot(sim)
end=#

###  Commands running from Julia REPL
numLoops, CPD, junctions, loops, k, L, σA, componentPhaseDirection = open_file("circuits/ac-i")
symbolic_assign()  #readline() does not work well in RELP
edit_param()       #readline() does not work well in RELP
new_model, u0 = build(Φ₀/10)
u0 = solve_init(u0)
sol = tspan("0.0", "1000")
single_plot("R1","V")

# =#

#=## Run from terminal commands
while true
    println(" --- Enter 'help' if you need help --- ")
    input = readline()
    if (input == "~")
        break
    elseif startswith(lowercase(input), "help")         #Give instructions on how to use
        help()
    elseif startswith(lowercase(input), "open_file")
        open_file(input[11:end-1])  
    elseif startswith(lowercase(input), "build")
        build(input[7:end-1])
    elseif startswith(lowercase(input), "solve_init")
        solve_init()
    elseif startswith(lowercase(input), "tspan")
        ts = split(strip(input[7:end-1]),',')   
        tspan(strip(ts[1]), strip(ts[2]))
    elseif startswith(lowercase(input), "single_plot")
        s_plot = split(input[13:end-1], ',')
        single_plot(strip(s_plot[1]), strip(s_plot[2]))
    elseif startswith(lowercase(input), "edit_param")
        edit_param()
    elseif startswith(lowercase(input), "symbolic_assign")
        symbolic_assign()
    end
end
# =#