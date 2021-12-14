function netlist(name)                  #USE JLD2.jl
    nl = open("$(name)_nl.csv", "w")            #Create file containing TOPOLOGY PARAMETERS OF SQUID APF CIRCUIT
    cp = open("$(name)_cp.txt", "w")            #Create file containing COMPONENT PARAMETERS OF SQUID APF CIRCUIT PARAMETER
    lc = open("$(name)_lc.csv", "w")            #Create file containing LOOPS AND COMPONENTS CONTAINED IN EACH LOOP

    write(nl, "NetList 1")                      #Title for netlist

    loops = []                                  #Stores components in each loop as array
    componentLoopDict = Dict()                  #Dictionary with components as keys and loops as values
    componentParamDict = Dict()                 #Dictionary with components as keys and parameters as values
    squidLoops = []                             #Maybe uneccessary 
    mutualInd = []                              #Stores data regarding which loops are mutally coupled
    junctions = []                              #Stores the names of the junctions
    extern_flux = []                            #Stores loops which have external flux

    println("Enter the number of loops in the circuit:")
    numLoops = readline()
    numLoops = parse(Int8, numLoops)             #Store number of Loops as an int

    println("Enter the Junctions:\n(Enter '~' when all are listed)")
    while true
        input = readline()
        if (input == "~")           
            break
        end
        push!(junctions, input)
    end

    println("Enter any mutually coupled loops:\nE.g. if loop 1 and 2 are coupled with mutual inductance 5μA/Φ𝜊 enter\n'1,2;5'\n(Enter '~' when all are listed)")
    while true
        input = readline()
        if (input == "~")           
            break
        end
        currentMutual = (parse(Int8, input[1]), parse(Int8, input[3]))    #Create tuple ### THIS DOESNT WORK FOR DECIMALS
        currentMutual = (currentMutual, parse(Int8,input[5]))            #Add tuple to tuple and add value
        push!(mutualInd, currentMutual)
    end

    for i in 1:numLoops                         #Asks about circuit elements
        push!(loops, [])
        write(nl, ",Loop $(i-1)")               #Write each loop name as a header
        write(lc, "Loop $(i-1)")
        println("Enter all components in Loop $(i-1) one by one\n(Enter '~' when all components are listed)")
        while true
            input = readline()
            if (input == "~")
                write(lc, "\n")                 
                break
            end
            push!(loops[i], input)
            write(lc, ",$input")                #Write each component to the loops+components.csv file
        end
    end

    for i in 1:numLoops                    #Iterate through all loops to find unique circuit elements      
        jj_count = 0
        for j in 1:length(loops[i])             #Iterate components in current loop
            if (loops[i][j][1] == 'J')             #Count Josephson Junctions in current loop
                jj_count = jj_count + 1
            end
            componentLoopDict[loops[i][j]]=push!(get(componentLoopDict, loops[i][j], []), i-1) #Creates dict with unique circuit elements
        end
        if (jj_count > 1)                       #If current loop has more than one Josephson Junction it is a SQUID loop... is this important?
            push!(squidLoops, i-1)
        end
    end

    for comp in keys(componentLoopDict)         #Finds circiut component parameters
        if (comp[1] == 'R')
            println("What is the resistance of $comp?")
            input = readline()
            componentParamDict[comp]=push!(get(componentParamDict, comp, []), input)
        elseif (comp[1] == 'C')
            println("What is the capacitance of $comp?")
            input = readline()
            componentParamDict[comp]=push!(get(componentParamDict, comp, []), input)
        elseif (comp[1] == 'L')
            println("What is the inductance of $comp?")
            input = readline()
            componentParamDict[comp]=push!(get(componentParamDict, comp, []), input)
        elseif (comp[1] == 'J')                                      #critical current, shunt resistance, and shunt capacitance needed for [Io] [G] [C]
            #=println("What is the critical current of $comp?")
            input = readline()
            componentParamDict[comp]=push!(get(componentParamDict, comp, []), input)
            println("What is the shunt resistance of $comp?")
            input = readline()
            componentParamDict[comp]=push!(get(componentParamDict, comp, []), input)
            println("What is the Stewart-McCumber parameter for $comp?")
            input = readline()
            componentParamDict[comp]=push!(get(componentParamDict, comp, []), input)=#
            println("What is the inductance for $comp?")
            input = readline()
            componentParamDict[comp]=push!(get(componentParamDict, comp, []), input)
        end
    end

    for j in 1:numLoops                         #First iteration to find the self and mutal inductance of all loops
        write(nl, "\nLoop $(j-1),/")
        for i in 2:numLoops                     #Second iteration to go through loop 1 to i for each and every loop 
    #                                           ^--- suspect matrix symmetry - better performance with "for i in j:numLoops"
            write(nl, ",")
            #SELF COUPLING
            for n in loops[i]                   #Third iteration to go through the components in loop 1
                if ((n[1] == 'J') || (n[1] == 'L')) #Merge if statemets??
                    if (j-1 in get(componentLoopDict, n, -1))   #If component n is in the loop j
                        if (i == j)             #Positive/Negative <--- Needs to be checked
                            write(nl, "+")
                        else
                            write(nl, "-")
                        end
                        write(nl, "$n")         #Parse first to generate numerical value?
                    end
                end
            end
            #MUTUAL COUPLING
            for n in mutualInd
                if ((i != j) && (i-1 in n[1]) && (j-1 in n[1]))
                    write(nl, "-$(n[2])")
                end
            end
        end 
    end

    for i in 1:length(junctions)
        write(nl, "\nJunction $i")
        for j in 1:numLoops
            if junctions[i] in loops[j]
                write(nl, ",1")
            else
                write(nl, ",0")
            end
        end
    end

    println("Enter the loops where external flux passes:\nE.g. if half the external flux passes through loop 1 and the remaining flux passes through loop 2 enter \n'1,0.5'\n'2,0.5'\n(Enter '~' when all are listed)")
    while true
        input = readline()
        if (input == "~")           
            break
        end
        input = split(input, ",")
        flux = (parse(Int8, input[1]), parse(Float16, input[2]))    #Create tuple only works for single digit int
        push!(extern_flux, flux)
    end

    write(nl, "\nFlux Φ,/")
    for i in 2:numLoops                                 #NEEDS WORK if external flux passes through two loops this does not work
        for j in 1:length(extern_flux)
            if (i-1 == extern_flux[j][1])
                write(nl, ",$(extern_flux[j][2])")
            else
                write(nl, ",0")
            end
        end
    end

    for (key, value) in componentParamDict
        write(cp, "$key: $value\n")
    end

    close(nl)
    close(cp)
    close(lc)
end

netlist("test")