using JLD2, FileIO, ModelingToolkit, LinearAlgebra

file = jldopen("test.jld2", "r")
PD = load("test.jld2", "editing/componentParamDict")
k = load("test.jld2", "matrices/k") #Needs to be 2 element vector
σB = load("test.jld2", "matrices/σB")
σA = load("test.jld2", "matrices/σA")
L = load("test.jld2", "matrices/L")
close(file)

k
σB
σA
L

Io = [0; 10e-6; 10e-6]
G = [1/27; 1/6.3; 1/6.3]
C = [0; 0.35; 0.35]
Φo = 2.06783383e-15
ΦE = 0.3*Φo

@variables t #=θ[1:3](t)=# I[1:3](t) #im[1:3](t) #In[1:3](t) # independent and dependent variables
im = map(first, [@variables $sym(t) for sym in (Symbol(:im, i) for i in 1:3)])
θ = map(first, [@variables $sym(t) for sym in (Symbol(:θ, i) for i in 1:3)])

D = Differential(t) # define an operator for the differentiation w.r.t. time
#eqs = Equation[]
eqs = [im[1] ~ 0]

for i in 1:2    #num of loops
    L_im = 0
    σA_θ = 0
    for j in 1:3    #num of juncs
        L_im = L_im + L[i,1:3][j]*im[j]
        σA_θ = σA_θ + σA[i,1:3][j]*θ[j]
    end
    L_im
    k_Φ = k[i]*ΦE
    push!(eqs, σA_θ ~ k_Φ - L_im)
end

for i in 1:3    #num of juncs
    σB_im = 0
    for j in 1:3    #num of loops
        σB_im = σB_im + σB[i,1:3][j]*im[j]
    end
    IoSθ = Io[i]*sin(θ[i])
    Gdθ = Φo/(2*pi)*G[i]*D(θ[i])
    Cdθ = Φo/(2*pi)*C[i]*D(D(θ[i]))
    push!(eqs, Gdθ + Cdθ ~ I[i] - IoSθ)
    push!(eqs, I[i] ~ σB_im)
end

eqs
eqs[1]
eqs[2]
eqs[3]
eqs[4]
eqs[5]
eqs[6]
eqs[7]
eqs[8]
eqs[9]


@named squid = ODESystem(eqs, t)

simple_squid = structural_simplify(squid)

simple_squid.eqs[1]
simple_squid.eqs[2]
simple_squid.eqs[3]
simple_squid.eqs[4]

u0 = [im[1] => 0.0, im[2] => 0.0, im[3] => 0.0, I[1] => 0.0, I[2] => 0.0, I[3] => 0.0, θ[1] => 0.0, θ[2] => 0.0, θ[3] => 0.0]

prob = ODEProblem(simple_squid, u0, (0.0,10.0))