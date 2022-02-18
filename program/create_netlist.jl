using JLD2, FileIO

function process_netlist(name)
    file = jldopen("$name.jld2", "a+")

    numLoops = read(file, "editing/numLoops")
    componentLoopDict = read(file, "editing/componentLoopDict")
    componentParamDict = read(file, "editing/componentParamDict")
    mutualInd = read(file, "editing/mutualInd")
    junctions = read(file, "editing/junctions")
    loops = read(file, "editing/loops")
    k = read(file, "matrices/k")

    close(file)

    file = jldopen("$name.jld2", "w") #Open file to write new data


    L = zeros(Float64, 0, numLoops)            
    σB = zeros(Float64, 0, numLoops)           
    componentPhaseDirection = Dict() 

    symbolDict = Dict()

    ### If symbolic ask for inductances in numerical values
    for comp in componentParamDict
        if (comp[1][1] == 'J')
            try 
                1+comp[2][4]
            catch e
                symbolDict[comp[2][4]] = push!(get(symbolDict, comp[2][4], []), comp[1])
            end
        elseif (comp[1][1] == 'L')
            try 
                comp[2]
                1+comp[2]
            catch e
                symbolDict[comp[2]] = push!(get(symbolDict, comp[2], []), comp[1])
            end
        end
    end
    for sym in symbolDict
        println("Please enter a value for $(sym[1])")
        input = readline()
        for i in sym[2]
            if (i[1] == 'J')
                componentParamDict[i][4] = parse(Float64, input)
            else
                componentParamDict[i] = parse(Float64, input)
            end
        end
    end

    ### Algorithm for finding L
    for j in 1:numLoops                         #First iteration to find the self and mutal inductance of all loops
        current_row = []
        for i in 2:numLoops                     #Second iteration to go through loop 1 to i for each and every loop 
    #                                           ^--- suspect matrix symmetry - better performance with "for i in j:numLoops"
            temp_float = 0.0
            #SELF COUPLING
            for n in loops[i]                   #Third iteration to go through the components in loop 1
                if ((n[1] == 'J') || (n[1] == 'L')) #Merge if statemets??
                    if (j-1 in get(componentLoopDict, n, -1))   #If component n is in the loop j

                        if (n[1] == 'J')
                            param = get(componentParamDict, n, 0)[4]
                        elseif (n[1] == 'L')
                            param = get(componentParamDict, n, 0)
                        end

                        if (i == j)             #Positive/Negative <--- Needs to be checked
                            temp_float = temp_float + parse(Float64, get(componentParamDict, n, 0))
                        else
                            temp_float = temp_float - parse(Float64, get(componentParamDict, n, 0))
                        end
                    end
                end
            end
            #MUTUAL COUPLING
            for n in mutualInd
                if ((i != j) && (i-1 in n[1]) && (j-1 in n[1]))
                    temp_float = temp_float - n[2]
                end
            end
            push!(current_row, temp_float)
        end 
        L = [L; current_row']
    end

    ### Algorithm for finding σB & σA & componentPhaseDirection
    for i in 1:length(junctions)
        current_row = []
        for j in 1:numLoops
            if junctions[i] in loops[j]
                junc_loops = get(componentLoopDict, junctions[i], -1)
                loop_count = 0
                for n in 1:length(junc_loops)
                    loop_count = loop_count + junc_loops[n]
                end
                if (loop_count/length(junc_loops) >= (j-1))
                    push!(current_row, 1)
                else
                    push!(current_row, -1)
                end
            else
                push!(current_row, 0)
            end
        end
        componentPhaseDirection[junctions[i]] = current_row
        σB = [σB; current_row']
    end

    L = transpose(L)
    σA = transpose(σB)


    #=
    println(L)
    println(σA)
    println(σB)
    println(componentPhaseDirection)
    =#

    file["editing/loops"] = loops
    file["editing/componentParamDict"] = componentParamDict
    file["editing/componentLoopDict"] = componentLoopDict
    file["editing/junctions"] = junctions
    file["editing/numLoops"] = numLoops
    file["editing/mutualInd"] = mutualInd
    file["matrices/k"] = k
    file["matrices/L"] = L
    file["matrices/σA"] = σA
    file["matrices/σB"] = σB
    file["matrices/componentPhaseDirection"] = componentPhaseDirection

    close(file)
end

function find_components(numLoops, loops)
    componentLoopDict = Dict()                  #Dictionary with components as keys and loops as values (used to find unique elements)
    componentParamDict = Dict()                 #Dictionary with components as keys and parameters as values 🟢
    junctions = []                              #Stores the names of the junctions

    for i in 1:numLoops                    #Iterate through all loops to find unique circuit components      
        for j in 1:length(loops[i])             #Iterate components in current loop
            componentLoopDict[loops[i][j]]=push!(get(componentLoopDict, loops[i][j], []), i-1) #Creates dict with unique circuit elements
        end
    end

    for comp in keys(componentLoopDict)         #Finds circiut component parameters
        if (comp[1] == 'V')
            push!(junctions, comp)


            println("Is $comp an AC or DC voltage source?")
            input = readline()
            if (lowercase(input) == "dc")
                println("What is the voltage of $comp (V)?")
                try
                    input = Meta.parse(readline())
                    componentParamDict[comp]=input
                catch e
                    if isa(e, MethodError)
                        input = readline()
                        componentParamDict[comp]=parse(Float64, input)
                    end
                end
            elseif (lowercase(input) == "ac")
                println("What is the amplitude and frequency of $comp (V, Hz)?")
                try
                    input = split(readline(), ',')
                    componentParamDict[comp]=(Meta.parse(input[1]), Meta.parse(input[2]))
                catch e
                    if isa(e, MethodError)
                        input = readline()
                        componentParamDict[comp]=(parse(Float64, input[1]), parse(Float64, input[2]))
                    end


                end
            end
        elseif (comp[1] == 'I')
            println("Is $comp an AC or DC current source?")
            input = readline()
            if (lowercase(input) == "dc")
                println("What is the current of $comp (A)?")
                try
                    input = Meta.parse(readline())
                    componentParamDict[comp]=input
                catch e
                    if isa(e, MethodError)
                        input = readline()
                        componentParamDict[comp]=parse(Float64, input)
                    end
                end
            elseif (lowercase(input) == "ac")
                println("What is the amplitude and frequency of $comp (A, Hz)?")
                try
                    input = split(readline(), ',')
                    componentParamDict[comp]=(Meta.parse(input[1]), Meta.parse(input[2]))
                catch e
                    if isa(e, MethodError)
                        input = readline()
                        componentParamDict[comp]=(parse(Float64, input[1]), parse(Float64, input[2]))
                    end
                end
            end
        elseif (comp[1] == 'R')
            push!(junctions, comp)
            println("What is the resistance of $comp (Ω)?")
            try
                input = Meta.parse(readline())
                componentParamDict[comp]=input
            catch e
                if isa(e, MethodError)
                    input = readline()
                    componentParamDict[comp]=parse(Float64, input)
                end
            end
        elseif (comp[1] == 'C')
            push!(junctions, comp)
            println("What is the capacitance of $comp (C)?")
            try
                input = Meta.parse(readline())
                componentParamDict[comp]=input
            catch e
                if isa(e, MethodError)
                    input = readline()
                    componentParamDict[comp]=parse(Float64, input)
                end
            end
        elseif (comp[1] == 'L')
            println("What is the inductance of $comp (H)?")
            try
                input = Meta.parse(readline())
                componentParamDict[comp]=input
            catch e
                if isa(e, MethodError)
                    input = readline()
                    componentParamDict[comp]=parse(Float64, input)
                end
            end
        elseif (comp[1] == 'J')
            push!(junctions, comp)

            println("What is the critical current of $comp (A)?")
            try
                input = Meta.parse(readline())
                componentParamDict[comp]=Any[input]
            catch e
                if isa(e, MethodError)
                    input = readline()
                    componentParamDict[comp]=Any[parse(Float64, input)]
                end
            end
            println("What is the shunt resistance of $comp (R)?")
            try
                input = Meta.parse(readline())
                componentParamDict[comp]=push!(get(componentParamDict, comp, []), input)
            catch e
                if isa(e, MethodError)
                    input = readline()
                    componentParamDict[comp]=push!(get(componentParamDict, comp, []), parse(Float64, input))
                end
            end
            println("What is the capacitance of $comp (F)?")
            try
                input = Meta.parse(readline())
                componentParamDict[comp]=push!(get(componentParamDict, comp, []), input)
            catch e
                if isa(e, MethodError)
                    input = readline()
                    componentParamDict[comp]=push!(get(componentParamDict, comp, []), parse(Float64, input))
                end
            end
            println("What is the inductance of $comp (H)?")
            try
                input = Meta.parse(readline())
                componentParamDict[comp]=push!(get(componentParamDict, comp, []), input)
            catch e
                if isa(e, MethodError)
                    input = readline()
                    componentParamDict[comp]=push!(get(componentParamDict, comp, []), parse(Float64, input))
                end
            end
        end
    end
    return componentLoopDict, componentParamDict, junctions
end

function new_netlist(name)
    file = jldopen("$name.jld2", "w")

    loops = []                                  #Stores components in each loop as array of array (MATRIX) 🟢
    #squidLoops = []                             #Maybe uneccessary 🔴
    mutualInd = []                              #Stores data regarding which loops are mutally coupled 🟢
    extern_flux = []                            #Stores loops which have external flux 🟢

    println("Enter the number of loops in the circuit:")
    numLoops = readline()
    numLoops = parse(Int8, numLoops)             #Store number of Loops as an int


    for i in 1:numLoops                         #Asks about circuit elements
        push!(loops, [])
        println("Enter all components in Loop $(i-1) one by one\n(Enter '~' when all components are listed)")
        while true
            input = readline()
            if (input == "~")  
                break
            end
            push!(loops[i], input)
        end
    end

    componentLoopDict, componentParamDict, junctions = find_components(numLoops, loops)

    #=println("Enter the Junctions:\n(Enter '~' when all are listed)\n --- This may not be necessary as im not sure if order of junctions matters ---")
    while true  ### ADD CHECK TO ENSURE USER INPUT JUNCTIONS EXIST IN THE CIRCUIT
        input = readline()
        if (input == "~")           
            break
        end
        push!(junctions, input)
    end=#
    

    println("Enter the external flux through each loop:\nE.g. if there are 3 loops (0, 1, 2) and 0.6 of the external flux passes through loop 1 and the remaining flux passes through loop 2 enter \n'0,0.6,0.4'")
    input = readline()
    if (input != "")
        flux = split(input, ',')
        for i in 1:length(flux)
            f = strip(flux[i])
            println(f)
            f = split(f, '/')
            if (length(f) == 1)
                try
                    push!(extern_flux, parse(Float64, f[1]))
                catch e
                    if isa(e, ArgumentError)
                        println(" --- Values must be passed as floats ---")
                    end
                end
            elseif (length(f) == 2)
                f1 = -1.0
                f2 = -1.0
                try
                    f1 = parse(Float64, f[1])
                    f2 = parse(Float64, f[2]) 
                catch e
                    if isa(e, ArgumentError)
                        println(" --- Values must be passed as floats ---")
                    end
                end
                if (f2 == 0)
                    println(" --- Cannot divide by 0 ---")
                else
                    push!(extern_flux, f1/f2)
                end
            end
        end
    end

    println("Enter any mutually coupled loops:\nE.g. if loop 1 and 2 are coupled with mutual inductance 5μA/Φ𝜊 enter\n'1,2,5'\n(Enter '~' when all are listed)")
    while true
        input = readline()
        if (input == "~")           
            break
        end
        try
            currentMutual = split(input, ',')
            if (length(currentMutual) != 3)
                throw(error)
            end
            mutualTuple = (parse(Int8, currentMutual[1]), parse(Int8, currentMutual[2]))
            mutualTuple = (mutualTuple, parse(Float64, currentMutual[3]))
            push!(mutualInd, mutualTuple)
        catch e
            if isa(e, ArgumentError)
                println(" --- One or more variable types are incorrect, loops must be Int and inductance must be Float ---")
            elseif e == error
                println(" --- Incorrect input length, enter 3 values seperated by commas ---")
            end
        end
    end

    file["editing/loops"] = loops
    file["editing/componentParamDict"] = componentParamDict
    file["editing/componentLoopDict"] = componentLoopDict
    file["editing/junctions"] = junctions
    file["editing/numLoops"] = numLoops
    file["editing/mutualInd"] = mutualInd
    file["matrices/k"] = extern_flux

    close(file)

    process_netlist(name)
end

function edit_netlist(name)
    file = jldopen("$name.jld2", "r+") #Open file to retrieve data

    numLoops = read(file, "editing/numLoops")
    componentLoopDict = read(file, "editing/componentLoopDict")
    componentParamDict = read(file, "editing/componentParamDict")
    mutualInd = read(file, "editing/mutualInd")
    junctions = read(file, "editing/junctions")
    loops = read(file, "editing/loops")
    k = read(file, "matrices/k")
    L = read(file, "matrices/L")
    σA = read(file, "matrices/σA")
    σB = read(file, "matrices/σB")
    componentPhaseDirection = read(file, "matrices/componentPhaseDirection")

    close(file)

    file = jldopen("$name.jld2", "w") #Open file to write new data
    
    while true
        println("What would you like to edit?")
        println(" --- Enter P to change component parameters ---")
        println(" --- Enter L to change components in loops ---")
        println(" --- Enter M to change mutally coupled loops ---")
        println(" --- Enter K to change external flux through loops ---")
        println(" --- Enter ~ when finished editing ---")
        input = readline()
        if (input == "~")
            break
        end
        if (uppercase(input) == "P")
            while true
                display(componentParamDict)
                println()
                println("Enter a component followed by '>' followed by it's new parameter. E.g. to change Ra from 2 Ohm to 3 Ohm enter Ra>3")
                println("For JJ's enter comma seperated parameters e.g. J1>10e-6,27,3.24,50e-3")
                println("Enter ~ when finished editing parameters")
                input = readline()
                if (input == "~")
                    break
                end
                comp_param = split(input, '>')
                if (comp_param[1] in keys(componentParamDict))
                    if startswith(comp_param[1], 'J')
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
                        componentParamDict[comp_param[1]] = param_a
                    elseif (startswith(comp_param[1], 'I') || startswith(comp_param[1], 'V'))
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
                            componentParamDict[comp_param[1]] = param_a
                        else
                            try
                                input = Meta.parse(comp_param[2])
                            catch e
                                if isa(e, MethodError)
                                    input = parse(Float64, input)
                                end
                            end
                            componentParamDict[comp_param[1]] = input
                        end
                    else
                        try
                            input = Meta.parse(comp_param[2])
                        catch e
                            if isa(e, MethodError)
                                input = parse(Float64, input)
                            end
                        end
                        componentParamDict[comp_param[1]] = input
                    end
                else
                    println("Component does not exist, try agian")
                end
            end
        elseif (uppercase(input) == "L")
            while true
                for i in 1:length(loops)
                    println("Loop $(i-1): $(loops[i])")
                end
                println("Which loop would you like to edit?\nEnter ~ when finished editing loops")
                input = readline()
                if (input == "~")
                    break
                end
                loop_num = parse(Int8, input)+1
                while true
                    println(loops[loop_num])
                    println("Enter a component in the loop to remove it, enter a new component to add it to the loop\nEnter ~ when finished editing loop $(loop_num-1)")
                    input = readline()
                    if (input == "~")
                        break
                    elseif (input in loops[loop_num])
                        deleteat!(loops[loop_num], findall(x->x==input,loops[loop_num]))
                    else
                        push!(loops[loop_num], input)
                    end
                end
            end
            componentLoopDict, componentParamDict, junctions = find_components(numLoops, loops)
        elseif (uppercase(input) == "M")
            while true
                display(mutualInd)
                println()
                println("Enter existing mutually coupled loops to remove, enter new mutually coupled loops to add\nEnter ~ when finished editing mutally coupled loops")
                input = readline()
                if (input == "~")
                    break
                end
                currentMutual = split(input, ',')
                mutualTuple = (parse(Int8, currentMutual[1]), parse(Int8, currentMutual[2]))
                mutualTuple = (mutualTuple, parse(Float64, currentMutual[3]))
                if (mutualTuple in mutualInd)
                    deleteat!(mutualInd, findall(x->x==mutualTuple,mutualInd))
                else
                    push!(mutualInd, mutualTuple)
                end
            end
        elseif (uppercase(input) == "K")
            display(k)
            println()

            println("Enter the new external flux through each loop:\nE.g. if there are 3 loops (0, 1, 2) and 0.6 of the external flux passes through loop 1 and the remaining flux passes through loop 2 enter \n'0,0.6,0.4'")

            input = readline()
            k = []
            flux = split(input, ',')
            for i in 1:length(flux)
                push!(k, parse(Float64, flux[i]))
            end
            display(k)
            println()
        end
    end

    file["editing/loops"] = loops
    file["editing/componentParamDict"] = componentParamDict
    file["editing/componentLoopDict"] = componentLoopDict
    file["editing/junctions"] = junctions
    file["editing/numLoops"] = numLoops
    file["editing/mutualInd"] = mutualInd
    file["matrices/k"] = k
    file["matrices/L"] = L
    file["matrices/σA"] = σA
    file["matrices/σB"] = σB
    file["matrices/componentPhaseDirection"] = componentPhaseDirection

    close(file)

    process_netlist(name)
end

println(" --- Enter 'E' to edit an existing netlist  --- ")
println(" --- Enter 'N' to create new netlist  --- ")
input = readline()
if (uppercase(input) == "E")
    println(" --- Enter filename (excluding '.jld2')  --- ")
    input = readline()
    edit_netlist("$input")
elseif (uppercase(input) == "N")
    println(" --- Enter filename (excluding '.jld2')  --- ")
    input = readline()
    new_netlist("$input")
end