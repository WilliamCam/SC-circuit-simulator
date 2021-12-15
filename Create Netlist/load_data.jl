using JLD2, FileIO
file = jldopen("test.jld2", "r")
display(file)

k = load("test.jld2", "matrices/k")
σB = load("test.jld2", "matrices/σB")
σA = load("test.jld2", "matrices/σA")
L = load("test.jld2", "matrices/L")

println("\n\nL:")
display(L)
println("\n\nσA:")
display(σA)
println("\n\nσB:")
display(σB)
println("\n\nk:")
display(k)
