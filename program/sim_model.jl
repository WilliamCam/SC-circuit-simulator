using JLD2, FileIO, ModelingToolkit, Plots, DifferentialEquations, LinearAlgebra, Statistics, Dates, Symbolics
include("model_builder.jl")
include("create_netlist.jl")

#Display function uses
function help()
    println()
    println("Use 'help()' to display this")
    println()
    println("Use 'open_file(name)' to open an existing .jdl2 file to extract circuit data")
    println()
    println("Use 'symbolic_assign()' to input numerical values for variables that are still in symbolic form")
    println("   - Symbolically assigned variables are temporary only, upon using net_edit() to change parameters symbolic variables return")
    println()
    println("Use 'net_edit(name)' to edit a netlist")
    println("   - After using net_edit(), open_file() must be called again to load correct data")
    println("   - Component parameters can also be edited by entering 'CPD[component] = new_parameter', (THIS IS TEMPORARY)")
    println()
    println("Use 'build(Φext)' to build the circuit based on the open file where Φext is the total external flux through the circuit")
    println()
    println("Use 'solve_init(old_u0, t_end)' to solve initial conditions for the model")
    println()
    println("Use 'sim_solve(u0, t_init, t_end)' to solve for a given timespan")
    println("   - A specific solver algorithm can be used for sim_solve() by adding solver=ROS3P() as a parameter for example")
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
    global numLoops = read(file, "editing/numLoops")
    global componentLoopDict = read(file, "editing/componentLoopDict")
    global CPD = read(file, "editing/componentParamDict")
    global mutualInd = read(file, "editing/mutualInd")
    global junctions = read(file, "editing/junctions")
    global loops = read(file, "editing/loops")
    global σA = read(file, "matrices/σA")
    #σB = read(file, "matrices/σB")
    global componentPhaseDirection = read(file, "matrices/componentPhaseDirection")
    close(file)
    print("File loaded \n circuit parameters:")
    for param in keys(CPD)
        print(string(param) * "\n")
    end
end

#Generates loop inductance matrix
function build_L_matrix()
    L = zeros(Num, numLoops,numLoops)
    ### Algorithm for finding L
    for j in 1:numLoops                                         #Iterate through all loops
        current_row = []
        for i in 1:numLoops                                     #Second iteration through all loops
            Lij = 0                                           #Float storing the value of the (j,i) position in matrix L
            #SELF COUPLING
            for n in loops[i]                                   #Iterate through components in loop i
                if ((n[1] == 'J') || (n[1] == 'L'))
                    if (j-1 in get(componentLoopDict, n, -1))   #If component n is also in the loop j
                        if (n[1] == 'J')
                            param = eval(Meta.parse(n * ".sys.L"))    #JJ case for setting param
                        elseif (n[1] == 'L')
                            eval(Meta.parse("@named " *n* " = build_inductance()"))    
                            param = eval(Meta.parse(n*".sys.L"))     #Inductor case for setting param
                        end
                        if (i == j)
                            Lij = Lij + param     #Adjust Lij by the value of the inductance of component n
                        else
                            Lij = Lij - param     #Adjust Lij by the value of the inductance of component n
                        end
                    end
                end
            end
            #MUTUAL COUPLING
            for n in mutualInd
                eval(Meta.parse("@named M" * string(n[1])*string(n[2]) * "= build_inductance()"))
                param = eval(Meta.parse("M" * string(n[1])*string(n[2]) * ".sys.L"))
                if ((i != j) && (i-1 in n[1]) && (j-1 in n[1])) #If the two currently observed loops are not the same loop and are stated as having mutual inductance
                    Lij = Lij - param              #Adjust Lij by the value of the mutual inductance
                end
            end
            push!(current_row, Lij)                      #Lij is pushed to current_row 
        end 
        #L = [L; current_row']
        L[j,:] = current_row'                                #current_row is pushed to the L matrix
    end

    return L
end

#Build the circuit based on the open file
function build_circuit()

    eqs = Equation[]                                        #Array to store equations

    built_loops = []
                                           #Array to store loops that have been built using MTK


    for i in 1:numLoops                                     #Iterate through all loops
        println("loop $(i)")                              #Display loop name
        current_loop = false
        c_source = ""
        for comp in loops[i]                                #Determine if a loop is a curernt source loop
            if startswith(comp, 'I')
                current_loop = true
                c_source = comp
            end
        end
        if (current_loop)                                   #Build a current source loop
            new_l = "@named loop$(i) = build_current_source_loop()"
            new_l = Meta.parse(new_l)
            new_l = eval(new_l)
        else                                                #Build a normal loop (no current source)
            new_l = "@named loop$(i) = build_loop()"
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
            new_c = "@named $j = build_resistor()"
            new_c = Meta.parse(new_c)                       #Using metaprogramming to generate components with unique names and parameters
            new_c = eval(new_c)
            built_components[j] = new_c                     #Push component to built components dictionary
        elseif (j[1] == 'C')                                #Built capacitor case
            param = get(CPD, j, 0)
            new_c = "@named $j = build_capacitor()"
            new_c = Meta.parse(new_c)
            new_c = eval(new_c)
            built_components[j] = new_c
        elseif (j[1] == 'V')                                #Built voltage source case
            param = get(CPD, j, 0)
            new_c = "@named $j = build_voltage_source()"
            new_c = Meta.parse(new_c)
            new_c = eval(new_c)
            built_components[j] = new_c
        elseif (j[1] == 'J')                                #Built Josephson Junction case
            params = get(CPD, j, 0)
            new_c = "@named $j = build_JJ()"
            new_c = Meta.parse(new_c)
            new_c = eval(new_c)
            built_components[j] = new_c
        end   
    end
    L = build_L_matrix()                                    #Generate inductance matrix

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
    
    add_loops!(eqs, built_loops, σA, built_components, L)
    current_flow(eqs, componentPhaseDirection, built_loops, built_components)
    
    
    

    @named _model  =  ODESystem(eqs, t)                     #Create an ODESystem with the existing equations

    @named model = compose(_model, sys)                     #Compose the existing ODESystem with the system states vector
    println()
    display(equations(model))                               #Display model equations
    println()
    display(states(model))                                  #Display model states
    println()
    display(parameters(model))
    println()
    new_model = structural_simplify(model);                 #structural_simplify Algorithm to improve performance
    return new_model, u0, built_components                               #Return structuraly simplified model and initial conditions
end

#Solve initial conditions
function solve_init(model, old_u0, t_end, p_string)
    ps = eval(Meta.parse(p_string))
    tspan_ini = (0.0, t_end)                                #Create a timespan
    prob = ODEProblem(model, old_u0, tspan_ini, ps, save_everystep = false, progress=true)  #Create an ODEProblem to solve for a specified time only saving the final component variable values
    sol = solve(prob, ROS3P())                                       #Solve the ODEProblem
    new_u0 = sol[:,end]                                     #Set the new initial conditions to the 
    return new_u0                                           #return the new intial conditions
end

#Give a timespan for the simulation and solve
function sim_solve(model, u0, t_init, t_end, ps_string)
    ps = eval(Meta.parse(ps_string))
    tspan = (t_init, t_end) #Create a timespan
    if (t_init != 0.0)
        u0 = solve_init(model,u0, tspan[1], ps_string)                           #Find the initial conditions for the start time
    end
    tsaves = LinRange(tspan[1],tspan[2], 50000)                 #Create tsaves
    prob = ODEProblem(model, u0, tspan, ps, saveat=tsaves, progress=true)   #Create an ODEProblem to solve for a specified time
    sol = solve(prob, maxiters=1e9, abstol = 1e-6)
    return sol                                                  #Return the solved ODEProblem
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
            #p = plot(ts[20000:end], vs[20000:end], label="$comp.$param", ylims=:round)                #Plot the component voltage vs time
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
#net_edit("circuits/ishaan")             #readline() in RELP has some issue where the first line is not read
# numLoops, CPD, junctions, loops, k, L, σA, componentPhaseDirection = open_file("circuits/ishaan")
# symbolic_assign()                       #readline() in RELP has some issue where the first line is not read
# new_model, u0 = build(3Φ₀);
# u0 = solve_init(u0, 1e-7)
# sol = sim_solve(u0, "0", "1e-7", ROS3P())

# single_plot("C3","V")
# single_plot("J4","V")

# ts = sol[t];  
# vs1 = Φ₀/2pi*sol[D(C3.sys.θ)];
# vs2 =  Φ₀/2pi*sol[D(J4.sys.θ)];
# plot(ts, vs1, label="C3.V");
# plot!(ts, vs2, label="J4.V");
# xlabel!("Time");
# ylabel!("Voltage")
