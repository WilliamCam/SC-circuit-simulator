include("create_netlist.jl")

#User input to enter edit or create functions
while true  
    println(" --- Enter 'E' to edit an existing netlist  --- ")
    println(" --- Enter 'N' to create new netlist  --- ")
    println(" --- Enter '~' to exit program  --- ")
    input = readline()
    if (uppercase(input) == "E")
        println(" --- Enter filename (excluding '.jld2')  --- ")
        input = readline()
        edit_netlist("$input")
    elseif (uppercase(input) == "N")
        println(" --- Enter filename (excluding '.jld2')  --- ")
        input = readline()
        new_netlist("$input")
    elseif (input == "~")
        exit()
    end
end