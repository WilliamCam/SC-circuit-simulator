using JLD2, FileIO, ModelingToolkit, Plots, DifferentialEquations, LinearAlgebra
include("model_builder.jl")

while true
    fail = true
    f_name = ""
    while fail
        println(" --- Enter filename (excluding '.jld2')  --- ")
        input = readline()
        f_name = "$input.jld2"
        fail = false
        try
            file = jldopen(f_name, "r")
            close(file)
        catch e
            if isa(e, SystemError)
                println("File does not exist")
                fail = true
            end
        end
    end
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

    #L
    #L[2,:]


    #CPD["Ib"] = 19e-6
    #CPD
    #CPD["J1"] = [1.0, 10, 3.2910597840193506e-14, 1]
    #CPD["V1"] = 1.0e-6
    #phi0 = 2.067833848e-15
    #I₀₁ = 10.0e-6
    #Gj1 = 1/10
    #βc₁ = 0.1
    #Cj1 = βc₁*phi0/(2*pi)*Gj1^2/(I₀₁)
    #Ma = phi0/(5e-6)

    j_len = length(junctions)

    eqs = Equation[]

    built_loops = []    #build loops
    #@named loop0 = build_current_source_loop(I = CPD["Ib"])
    #push!(built_loops, loop0)
    external_flux_strength = 3*Φ₀/14

    for i in 1:numLoops
        println("loop $(i-1)")
        if ("Ib" in loops[i])
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

    #display(built_components)

    #=loops_built_components = []
    for i in 1:numLoops
        current_components = Component[]
        for j in loops[i]
            push!(current_components, get(built_components, j, 0))
        end
        push!(loops_built_components, current_components)
    end=#

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
    

    sys = Vector{ODESystem}(old_sys)

    #=eqs_t = []
    push!(eqs_t, built_loops[2].sys.Φₗ ~ dot(L[2,:], [l.sys.iₘ for l in built_loops]))
    eqs_t
    L
    L[:,2]
    built_loops[1].sys.iₘ
    built_loops[2].sys.iₘ=#

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

    tspan_ini = (0.0,1e-9)

    prob = ODEProblem(new_model, my_u0, tspan_ini, save_everystep = false, progress=true)
    sol = solve(prob, ROS3P())
    u0 = sol[:,end]

    println("\n --- Enter a timespan --- ")
    println(" --- e.g. '0.0, 100.0' or '1e-8, 1e-7' --- ")
    input = readline()
    tspan = (parse(Float64, split(input, ',')[1]), parse(Float64, split(input, ',')[2]))
    tsaves = LinRange(tspan[1],tspan[2], 50000)
    prob = ODEProblem(new_model, u0, tspan, saveat=tsaves, progress=true)
    sol = solve(prob)

    while true
        println(" --- Enter a component and its variable to plot, enter ~ to end ---")
        println(" --- e.g. to plot current through J1 enter J1,i or to plot voltage through J2 enter J2,V ---")
        input = readline()
        if (input == "~")
            break
        end
        try
            input = split(input, ',')
            c = input[1]
            v = input[2]
            if (v == "i")
                str = c*".sys.i"
                ex = Meta.parse(str)
                p = plot(sol, vars=[eval(ex)])
                png("$(f_name)_$str")
            elseif (v == "V")
                str = "D($c.sys.θ)"
                ex = Meta.parse(str)
                p = plot(sol, vars=[eval(ex)])
                png("$(f_name)_$str")
            end
        catch e
            if isa(e, UndefVarError)
                println(" --- Component does not exist ---")
            elseif isa(e, ArgumentError)
                println(" --- Cannot plot voltage through resistors at this point in time ---")
            end
        end
    end

    println(" --- Enter ~ to end program, otherwise press any key to start a new circuit simulation --- ")
    input = readline()
    if (input == "~")
        break
    end
end

#sim()

#=
extern_flux = []

numLoops = 3

println("Enter the external flux through each loop:\nE.g. if there are 3 loops (0, 1, 2) and 0.6 of the external flux passes through loop 1 and the remaining flux passes through loop 2 enter \n'0.6-0.4'\nDo NOT enter the 0th (Ib) Loop")
    input = "5"
    try
        flux = split(input, '-')
        if (length(flux) != numLoops-1)

        end
        for i in 1:length(flux)
            push!(extern_flux, parse(Float64, flux[i]))
        end
    catch e
        
    end
    =#