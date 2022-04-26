using DifferentialEquations, Plots, Statistics

include("C://Users//21958742//GitHub//SC-circuit-simulator//program//SCsim.jl")
scsim.open_file("C://Users//21958742//GitHub//SC-circuit-simulator//examples//Quantum Gravity//QuartzTnetwork.jld2")

model, u0, old_model = scsim.build_circuit()

ps = [
    scsim.V1.sys.V => 1.0
    scsim.V1.sys.ω => 5.36e+7
    scsim.R1.sys.R=> 50.0
    scsim.C1.sys.C => 10.0e-12
    scsim.Cq.sys.C => 8.7e-17
    scsim.Lq.sys.L => 4.0
    scsim.Lq.sys.k => 1.0
    scsim.Rq.sys.R => 500.0
    scsim.C0.sys.C=> 6.9e-11
    scsim.C2.sys.C => 10.0e-12
    scsim.R2.sys.R=>1.0e+3
]

tspan = (0.0, 1.0e-3)
tsaves = LinRange(tspan[2]/10.0, tspan[2], 10000)
dt = tsaves[2] - tsaves[1]

prob = ODEProblem(model, u0, tspan, ps, maxiters=1e9, force_dtmin=true, saveat=tsaves, abstol=1e-8)
@time sol = solve(prob, ROS3P())
v =  1/dt * scsim.Φ₀/(2*pi)*diff(sol[scsim.R2.sys.θ])
i = sol[scsim.Rq.sys.i][2:end]
plot(i)
fs = 2/dt
using FFTW
v_fft = abs.(2*rfft(v)/length(v))
fn = LinRange(0.0, fs/2, length(v_fft))
#plot(fn, v_fft)
maximum(v_fft)
function RMS(x)
    return sqrt(mean((x .- mean(x)).^2))
end
RMS(v)

ωλ = 1/sqrt(prob.p[12]*prob.p[2])

Ntrajectories = 1000

function prob_func(prob, i ,repeat)
    prob.p[7] = ω_vec[i]
    prob
end

function output_func(sol,i)
    push!(logger, 1)
    println(string(Ntrajectories-length(logger)))
    (RMS(sol[scsim.R2.sys.i]),false)
end

function RMS_IRq(sol,i)
    push!(logger, 1)
    println(string(Ntrajectories-length(logger)))
    ([RMS(sol[scsim.Rq.sys.i]), RMS(sol[scsim.R2.sys.i])],false)
end

##Ensemble sims##
ω_vec = ωλ .+ LinRange(-5000,5000, Ntrajectories)
ps[1] = scsim.V1.sys.V => 1.0 

prob = ODEProblem(model, u0, tspan, ps, saveat = tsaves, maxiters=1e9, force_dtmin=true, dense = false,
    abstol = 1e-8    
)
logger = []
ensemble_prob = EnsembleProblem(prob,prob_func=prob_func, output_func=RMS_IRq)
@time sim_0p9V = solve(ensemble_prob,Rodas5(),EnsembleThreads(),trajectories=Ntrajectories)

plot(ω_vec/(2*pi), sim_0p9V[2,:])

ps[1] = scsim.V1.sys.V => 1.2 

prob = ODEProblem(model, u0, tspan, ps, saveat = tsaves, maxiters=1e9, progress=true, progress_steps = 1, force_dtmin=true, dense = false)
ensemble_prob = EnsembleProblem(prob,prob_func=prob_func, output_func=RMS_IRq)
@time sim_1p2V = solve(ensemble_prob,Rodas5(),EnsembleSerial(),trajectories=Ntrajectories)

plot(ω_vec/(2*pi), sim_1p2V[2,:])