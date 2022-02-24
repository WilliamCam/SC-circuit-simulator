using JLD2, FileIO, ModelingToolkit, Plots, DifferentialEquations, LinearAlgebra, Statistics, Dates
include("model_builder.jl")

#Display function uses
function help()
    println()
    println("Use 'help()' to display this")
    println()
    println("Use 'open_file(name)' to open an existing .jdl2 file to extract circuit data")
    println()
    println("Use 'symbolic_assign()' to input numerical values for variables that are still in symbolic form")
    println()
    println("Use 'edit_param()' to edit Component parameters")
    println("   - Component parameters can also be edited by entering 'CPD[component] = new_parameter'")
    println("   WARNING: Changing a components inductance in sim_model.jl will not work as the inductance matrix (L) is not re-evaluated")
    println()
    println("Use 'build(Φext)' to build the circuit based on the open file where Φext is the total external flux through the circuit")
    println()
    println("Use 'solve_init(old_u0, t_end)' to solve initial conditions for the model")
    println()
    println("Use 'tspan(u0, t_init, t_end)' to solve for a given timespan")
    println()
    println("Use 'single_plot(comp, param)' to plot the voltage or current through a component")
    println("   - To plot voltage enter 'V' as the param")
    println("   - To plot current enter 'i' as the param")
end

#Open an existing .jdl2 file to extract circuit data
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

#Ask user to input numerical values for variables that are still in symbolic form
function symbolic_assign()  #Does not work for AC voltage/current sources
    symbolDict = Dict()
    for comp in CPD
        if (comp[1][1] == 'J')
            for i in 1:4
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

#Edit parameters while in program
function edit_param()
    while true
        display(CPD)         #Display the current CPD
        println()
        println("Enter a component followed by '>' followed by it's new parameter. E.g. to change Ra from 2 Ohm to 3 Ohm enter Ra>3")
        println("For JJ's enter comma seperated parameters e.g. J1>10e-6,27,3.24,50e-3")
        println("Enter ~ when finished editing parameters")
        input = readline()
        if (input == "~")
            break                           #Exit edit component parameters
        end
        comp_param = split(input, '>')
        if (comp_param[1] in keys(CPD))  #Check if component exists
            if startswith(comp_param[1], 'J')           #Josephson Junction case
                params = split(comp_param[2], ',')
                param_a = []
                for p in params
                    try
                        input = Meta.parse(p)
                    catch e
                        if isa(e, MethodError)
                            input = parse(Float64, p)
                        end
                    end
                    push!(param_a, input)
                end
                CPD[comp_param[1]] = param_a
            elseif (startswith(comp_param[1], 'I') || startswith(comp_param[1], 'V')) #Current/Voltage source case
                params = split(comp_param[2], ',')
                if (length(params) == 2)
                    param_a = []
                    for p in params
                        try
                            input = Meta.parse(p)
                        catch e
                            if isa(e, MethodError)
                                input = parse(Float64, p)
                            end
                        end
                        push!(param_a, input)
                    end
                    CPD[comp_param[1]] = param_a
                else
                    try
                        input = Meta.parse(comp_param[2])
                    catch e
                        if isa(e, MethodError)
                            input = parse(Float64, input)
                        end
                    end
                    CPD[comp_param[1]] = input
                end
            else                                        #All other components
                try
                    input = Meta.parse(comp_param[2])
                catch e
                    if isa(e, MethodError)
                        input = parse(Float64, input)
                    end
                end
                CPD[comp_param[1]] = input
            end
        else
            println("Component does not exist, try agian")
        end
    end
end

#Build the circuit based on the open file
function build(Φext)
    eqs = Equation[]                                        #Array to store equations

    built_loops = []                                        #Array to store loops that have been built using MTK

    for i in 1:numLoops                                     #Iterate through all loops
        println("loop $(i-1)")                              #Display loop name
        current_loop = false
        c_source = ""
        for comp in loops[i]                                #Determine if a loop is a curernt source loop
            if startswith(comp, 'I')
                current_loop = true
                c_source = comp
            end
        end
        if (current_loop)                                   #Build a current source loop
            param = CPD[c_source]
            if (length(param) == 1)                         #DC current source case
                new_l = "@named loop$(i-1) = build_current_source_loop(I = $(param))"
            else                                            #AC current source case
                new_l = "@named loop$(i-1) = build_current_source_loop(I = $(param[1]), ω = $(1/param[2]))"
            end
            new_l = Meta.parse(new_l)
            new_l = eval(new_l)
        else                                                #Build a normal loop (no current source)
            new_l = "@named loop$(i-1) = build_loop(Φₑ = $(Φext*k[i]))"
            new_l = Meta.parse(new_l)                       #Using metaprogramming to generate loops with unique names and parameters
            new_l = eval(new_l)
        end
        push!(built_loops, new_l)                           #Push built loop to built_loops array
    end

    built_components = Dict()                               #Dictionary to store components that have been built with MTK
    for j in junctions                                      #Iterate through juncrions
        println(j)                                          #Display junction name
        if (j[1] == 'R')                                    #Built resistor case
            param = get(CPD, j, 0)
            new_c = "@named $j = build_resistor(R = $param)"
            new_c = Meta.parse(new_c)                       #Using metaprogramming to generate components with unique names and parameters
            new_c = eval(new_c)
            built_components[j] = new_c                     #Push component to built components dictionary
        elseif (j[1] == 'C')                                #Built capacitor case
            param = get(CPD, j, 0)
            new_c = "@named $j = build_capacitor(C = $param)"
            new_c = Meta.parse(new_c)
            new_c = eval(new_c)
            built_components[j] = new_c
        elseif (j[1] == 'V')                                #Built voltage source case
            param = get(CPD, j, 0)
            if (length(param) == 1)                         #DC voltage case
                new_c = "@named $j = build_voltage_source(V = $param)"
            else                                            #AC voltage case
                new_c = "@named $j = build_voltage_source(V = $(param[1]), ω = $(1/param[2]))"
            end
            new_c = Meta.parse(new_c)
            new_c = eval(new_c)
            built_components[j] = new_c
        elseif (j[1] == 'J')                                #Built Josephson Junction case
            params = get(CPD, j, 0)
            new_c = "@named $j = build_JJ(I0 = $(params[1]), R = $(params[2]), C = $(params[3]))"
            new_c = Meta.parse(new_c)
            new_c = eval(new_c)
            built_components[j] = new_c
        end
    end

    old_sys = []                                            #Array to store system states                                          
    u0 = Pair{Num, Float64}[]                               #Array to store system initial condionts  (Set to 0)

    for comp in built_components                            #Iterate through components to find component system states and intial conditons
        push!(old_sys, comp[2].sys)
        if (comp[1][1] != 'V')
            push!(u0, comp[2].sys.θ=>0.0)                   #θ initialised to 0
            push!(u0, comp[2].sys.i=>0.0)                   #i initialised to 0
            if (uppercase(comp[1][1]) in ['C', 'J', 'V'])   
                push!(u0, D(comp[2].sys.θ)=>0.0)            #D(θ) initialised to 0 for capacitors, JJs and voltage sources
            end
        end
    end
    for loop in built_loops                                 #Iterate through components to find loop system states
        push!(old_sys, loop.sys)
    end
    sys = Vector{ODESystem}(old_sys)                        #Convert system states array to an ODESystem vector form

    #Functions from model_builder.jl to form appropriate equations
    inductance(eqs, L, built_loops)
    current_flow(eqs, componentPhaseDirection, built_loops, built_components)
    for i in 1:numLoops
        add_loop!(eqs, built_loops[i], σA[i,:], built_components)
    end

    @named _model  =  ODESystem(eqs, t)                     #Create an ODESystem with the existing equations

    @named model = compose(_model, sys)                     #Compose the existing ODESystem with the system states vector

    println()
    display(equations(model))                               #Display model equations
    println()
    display(states(model))                                  #Display model states
    println()
    new_model = structural_simplify(model)                  #structural_simplify Algorithm to improve performance
    return new_model, u0                                    #Return structuraly simplified model and initial conditions
end

#Solve initial conditions
function solve_init(old_u0, t_end)
    tspan_ini = (0.0, t_end)                                #Create a timespan
    prob = ODEProblem(new_model, old_u0, tspan_ini, save_everystep = false, progress=true)  #Create an ODEProblem to solve for a specified time only saving the final component variable values
    sol = solve(prob)                                       #Solve the ODEProblem
    new_u0 = sol[:,end]                                     #Set the new initial conditions to the 
    return new_u0                                           #return the new intial conditions
end

#Give a timespan for the simulation and solve
function tspan(u0, t_init, t_end)
    new_u0 = solve_init(u0, t_init)                         #Find the initial conditions for the start time
    tspan = (parse(Float64, t_init), parse(Float64, t_end)) #Create a timespan
    tsaves = LinRange(tspan[1],tspan[2], 50000)             #Create tsaves
    prob = ODEProblem(new_model, new_u0, tspan, saveat=tsaves, progress=true)   #Create an ODEProblem to solve for a specified time
    sol = solve(prob)                                       #Solve the ODEProblem
    return sol                                              #Return the solved ODEProblem
end

#Plot a single component
function single_plot(comp, param)
    try
        if (param == "i")                                   #Current on y axis
            str = "$comp.sys.i"                             #Using metaprogramming
            ex = Meta.parse(str)
            p = plot(sol, vars=[eval(ex)], label="$comp.$param", ylims=:round)  #Plot the component curernt vs time
            xlabel!("Time")                                 #Add axis labels
            ylabel!("Current")
            display(p)                                      #Display plot (if suppourted)
            png("$(now())")                                 #Save plot as the current date and time
        elseif (param == "V")                               #Voltage on y axis
            str = "D($comp.sys.θ)"                          #Using metaprogramming
            ex = Meta.parse(str)
            vs = Φ₀/2pi*sol[eval(ex)]                       #Multiply by constant (conversion from θ to V)
            ts = sol[t]                                     #Timespan
            p = plot(ts, vs, label="$comp.$param", ylims=:round)                #Plot the component voltage vs time
            xlabel!("Time")                                 #Add axis labels
            ylabel!("Voltage")
            display(p)                                      #Display plot (if suppourted)
            png("$(now())")                                 #Save plot as the current date and time
        end
    catch e
        if isa(e, UndefVarError)                            #If the component cannot be found print a suitable error message
            println(" --- Component does not exist ---")
        elseif isa(e, ArgumentError)                        #If attempt to plot voltage through a component where D(θ) is not saved
            dtime = sol.t[2]-sol.t[1]                       #Find change in time
            str = "$comp.sys.θ"                             #Using metaprogramming
            ex = Meta.parse(str)
            dv = sol[eval(ex)]
            res = []                                        #Result vector
            for i in 1:length(dv)-1                         #Differentiation by hand
                d = Φ₀/2pi*(dv[i+1] - dv[i])/dtime          #Change in θ over change in time multiplied by constant (conversion from θ to V)
                push!(res, d)                               #Push dθ/dt to result vector
            end
            p = plot(sol.t[1:end-1], res, label="$comp.$param", ylims=:round)   #Plot the component voltage vs time
            xlabel!("Time")                                 #Add axis labels
            ylabel!("Current")
            display(p)                                      #Display plot (if suppourted)
            png("$(now())")                                 #Save plot as the current date and time
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
help()
numLoops, CPD, junctions, loops, k, L, σA, componentPhaseDirection = open_file("ishaan")
symbolic_assign()  #readline() in RELP has some issue where the first line is not read
edit_param()       #readline() in RELP has some issue where the first line is not read
new_model, u0 = build(3Φ₀)
u0 = solve_init(u0, 1e-6)
sol = tspan(u0, "0.0", "1e-9")
single_plot("C3","V")
single_plot("J4","V")

#=## Run from terminal commands
println(" --- Enter 'help' if you need help --- ")
println(" --- Enter '~' to exit program  --- ")
while true
    input = readline()
    if (input == "~")
        exit()
    elseif startswith(lowercase(input), "help")         #Give instructions on how to use
        help()
    elseif startswith(lowercase(input), "open_file")
        open_file(input[11:end-1])  
    elseif startswith(lowercase(input), "build")        #Pass through data from open_file()?
        build(input[7:end-1])
    elseif startswith(lowercase(input), "solve_init")   #Pass through model form build()?
        si = split(strip(input[12:end-1]),',')
        solve_init(strip(si[1]), strip(si[2]))
    elseif startswith(lowercase(input), "tspan")        #Pass through model form build()?
        ts = split(strip(input[7:end-1]),',')   
        tspan(strip(ts[1]), strip(ts[2]), strip(ts[3]))
    elseif startswith(lowercase(input), "single_plot")  #Pass through solution form tspan()?
        s_plot = split(input[13:end-1], ',')
        single_plot(strip(s_plot[1]), strip(s_plot[2]))
    elseif startswith(lowercase(input), "edit_param")   #Pass through CPD?
        edit_param()
    elseif startswith(lowercase(input), "symbolic_assign")  #Pass through CPD?
        symbolic_assign()
    end
end
# =#