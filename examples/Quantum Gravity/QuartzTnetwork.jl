using DifferentialEquations, Plots, Statistics

include("C://Users//21958742//GitHub//SC-circuit-simulator//program//SCsim.jl")
scsim.open_file("C://Users//21958742//GitHub//SC-circuit-simulator//examples//Quantum Gravity//QuartzTnetwork-cs.jld2")

model, u0, old_model = scsim.build_circuit()

ps = [
    scsim.loop1.sys.I => 1e-4
    scsim.loop1.sys.ω => 5.360562674188974e7
    scsim.R1.sys.R=> 50.0
    scsim.C1.sys.C => 10.0e-12
    scsim.Cq.sys.C => 8.7e-17
    scsim.Lq.sys.L => 4.0
    scsim.Lq.sys.k => 1.0
    scsim.Rq.sys.R => 5.36
    scsim.C0.sys.C=> 6.9e-11
    scsim.C2.sys.C => 10.0e-12
    scsim.R2.sys.R=>50.0
]

tspan = (0.0, 5.0e-3)
tsaves = LinRange(tspan[2]/10.0, tspan[2], 10000)
dt = tsaves[2] - tsaves[1]
using ModelingToolkit
dae_model = dae_index_lowering(model)
prob = ODAEProblem(dae_model, u0, tspan, ps, saveat = tsaves, maxiters=1e9)
#prob = ODEProblem(model, u0, tspan, ps, saveat=tsaves, maxiters=1e9)
@time sol = solve(prob, Tsit5(), abstol = 1e-8)
v =  1/dt * scsim.Φ₀/(2*pi)*diff(sol[scsim.Rq.sys.i])
i = sol[scsim.Rq.sys.i][2:end]
plot(i)
fs = 2/dt
using FFTW
v_fft = abs.(2*rfft(v)/length(v))
fn = LinRange(0.0, fs/2, length(v_fft))
plot(fn, v_fft)
maximum(v_fft)
function RMS(x)
    return sqrt(mean((x .- mean(x)).^2))
end
RMS(v)


Cindex = findfirst(isequal(scsim.Cq.sys.C),parameters(model))
Lindex = findfirst(isequal(scsim.Lq.sys.L),parameters(model))
ωλ = 1/ sqrt(prob.p[9]*prob.p[2])
Ntrajectories = 1000

function prob_func(prob, i ,repeat)
    prob.p[21] = ω_vec[i]
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
ω_vec = ωλ .+ LinRange(-5e+3,5e+3, Ntrajectories)

ps[1] = scsim.loop1.sys.I => 0.5e-3 

logger = []
prob = ODAEProblem(model, u0, tspan, ps, saveat = tsaves, maxiters=1e9, force_dtmin=true, dense = false, abstol = 1e-8)
ensemble_prob = EnsembleProblem(prob,prob_func=prob_func, output_func=RMS_IRq)
@time sim_0p5mA = solve(ensemble_prob,Tsit5(),EnsembleThreads(),trajectories=Ntrajectories)

ps[1] = scsim.loop1.sys.I => 0.6e-3 

logger = []
prob = ODAEProblem(model, u0, tspan, ps, saveat = tsaves, maxiters=1e9, force_dtmin=true, dense = false, abstol = 1e-8)
ensemble_prob = EnsembleProblem(prob,prob_func=prob_func, output_func=RMS_IRq)
@time sim_0p6mA = solve(ensemble_prob,Tsit5(),EnsembleSerial(),trajectories=Ntrajectories)

ps[1] = scsim.loop1.sys.I => 0.7e-3 

logger = []
prob = ODAEProblem(model, u0, tspan, ps, saveat = tsaves, maxiters=1e9, force_dtmin=true, dense = false, abstol = 1e-8)
ensemble_prob = EnsembleProblem(prob,prob_func=prob_func, output_func=RMS_IRq)
@time sim_0p7mA = solve(ensemble_prob,Tsit5(),EnsembleSerial(),trajectories=Ntrajectories)

ps[1] = scsim.loop1.sys.I => 0.8e-3 

logger = []
prob = ODAEProblem(model, u0, tspan, ps, saveat = tsaves, maxiters=1e9, force_dtmin=true, dense = false, abstol = 1e-8)
ensemble_prob = EnsembleProblem(prob,prob_func=prob_func, output_func=RMS_IRq)
@time sim_0p8mA = solve(ensemble_prob,Tsit5(),EnsembleSerial(),trajectories=Ntrajectories)

plot(ω_vec/(2*pi),sim_0p5mA[2,:])

ps[1] = scsim.loop1.sys.I => 1.0e-3 

logger = []
prob = ODAEProblem(model, u0, tspan, ps, saveat = tsaves, maxiters=1e9, force_dtmin=true, dense = false, abstol = 1e-8)
ensemble_prob = EnsembleProblem(prob,prob_func=prob_func, output_func=RMS_IRq)
@time sim_1mA = solve(ensemble_prob,Tsit5(),EnsembleThreads(),trajectories=Ntrajectories)

ps[1] = scsim.loop1.sys.I => 1.25e-3 

logger = []
prob = ODAEProblem(model, u0, tspan, ps, saveat = tsaves, maxiters=1e9, force_dtmin=true, dense = false, abstol = 1e-8)
ensemble_prob = EnsembleProblem(prob,prob_func=prob_func, output_func=RMS_IRq)
@time sim_1p25mA = solve(ensemble_prob,Tsit5(),EnsembleThreads(),trajectories=Ntrajectories)

ps[1] = scsim.loop1.sys.I => 2.0e-3 

logger = []
prob = ODAEProblem(model, u0, tspan, ps, saveat = tsaves, maxiters=1e9, force_dtmin=true, dense = false, abstol = 1e-8)
ensemble_prob = EnsembleProblem(prob,prob_func=prob_func, output_func=RMS_IRq)
@time sim_2mA = solve(ensemble_prob,Tsit5(),EnsembleThreads(),trajectories=Ntrajectories)

ps[1] = scsim.loop1.sys.I => 3.0e-3 

logger = []
prob = ODAEProblem(model, u0, tspan, ps, saveat = tsaves, maxiters=1e9, force_dtmin=true, dense = false, abstol = 1e-8)
ensemble_prob = EnsembleProblem(prob,prob_func=prob_func, output_func=RMS_IRq)
@time sim_3mA = solve(ensemble_prob,Tsit5(),EnsembleThreads(),trajectories=Ntrajectories)

 plt1 = plot(ω_vec/(2*pi), 
    [sim_0p5mA[1,:], sim_0p6mA[1,:],sim_0p7mA[1,:], sim_0p8mA[1,:], sim_1mA[1,:], sim_1p25mA[1,:], sim_2mA[1,:], sim_3mA[1,:]],
    xlabel = "Frequency (Hz)",
    ylabel = "Motional current (Irms)",
    label = ["Iin = 1mA" "Iin = 1.5mA" "Iin = 2mA" "Iin = 3mA"],
    title = "Nonlinear motional inductor k =1"

)

plt2 = plot(ω_vec/(2*pi),
 50 .* [sim_0p5mA[2,:]/5.0, sim_0p6mA[2,:]/6.0,sim_0p7mA[2,:]/7.0, sim_0p8mA[2,:]/8.0,
    sim_1mA[2,:]/10.0, sim_1p25mA[2,:]/12.5, sim_2mA[2,:]/20.0, sim_3mA[2,:]/30.0],
xlabel = "Frequency (Hz)",
ylabel = "Output Voltage (Vrms)",
label = ["Iin = 1mA" "Iin = 1.5mA" "Iin = 2mA" "Iin = 3mA"],
title = "Nonlinear motional inductor k =1"

)

savefig(plt1, "C://Users//21958742//GitHub//SC-circuit-simulator//examples//Quantum Gravity//TnetworkSim-motionalI-k1.0-0p5mA.pdf")
savefig(plt2, "C://Users//21958742//GitHub//SC-circuit-simulator//examples//Quantum Gravity//TnetworkSim-Vout-k1.0-0p5mA.pdf")
## Harmonic Balance ##
using Symbolics
using HarmonicBalance

eqs = equations(model)
obs = observed(model)

eqs1 = substitute(eqs, Dict(scsim.C2.sys.k => 0.0, scsim.C1.sys.k => 0.0, scsim.C0.sys.k=>0.0,
    scsim.Cq.sys.k=>0.0, scsim.R1.sys.k=>0.0, scsim.R2.sys.k=>0.0, scsim.Rq.sys.k=>0.0))

sub2 = Symbolics.solve_for(obs[8], scsim.Cq.sys.i)
eqs2 = substitute(eqs1, Dict(scsim.Cq.sys.i => sub2))

sub3 = Symbolics.solve_for(obs[7], scsim.Rq.sys.i)
eqs3 = substitute(eqs2, Dict(scsim.Rq.sys.i => sub3))

sub4 = Symbolics.solve_for(obs[3], scsim.C0.sys.i)
eqs4 = substitute(eqs3, Dict(scsim.C0.sys.i => sub4))

sub5 = Symbolics.solve_for(obs[4], scsim.R1.sys.i)
eqs5 = substitute(eqs4, Dict(scsim.R1.sys.i => sub5))

##eliminate equation 13 by substitution
sub6 = Symbolics.solve_for(eqs5[13], scsim.R2.sys.i)
eqs6 = substitute(eqs5[1:12], Dict(scsim.R2.sys.i => sub6))


##eliminate equation 5 by subtitution
sub7 = -d(scsim.Cq.sys.θ, scsim.t, 2) * 2.067833848e-15*scsim.Cq.sys.C / 6.283185307179586 + (scsim.loop1.sys.I*cos(scsim.loop1.sys.ω*scsim.t))
eqs7 = substitute([eqs6[1:4]; eqs6[6:end]], Dict(scsim.loop2.sys.iₘ => sub7))

diff_eqs = [
        d(scsim.C1.sys.θ, scsim.t, 2) ~ eqs7[1].rhs
        d(scsim.C2.sys.θ, scsim.t, 2) ~ (6.283185307179586*((3.2910597840193497e-16*d(scsim.C2.sys.θ,scsim.t) - 3.2910597840193497e-16*d(scsim.C0.sys.θ,scsim.t)) / (-scsim.R2.sys.R))) / (2.067833848e-15*scsim.C2.sys.C)
        eqs7[6]
        d(scsim.C0.sys.θ, scsim.t, 2) ~ (6.283185307179586*(scsim.loop1.sys.I*cos(scsim.loop1.sys.ω*scsim.t) + (3.2910597840193497e-16*d(scsim.C0.sys.θ, scsim.t) - 3.2910597840193497e-16*d(scsim.C2.sys.θ, scsim.t)) / (-scsim.R2.sys.R) - 3.29105978401935e-16*scsim.Cq.sys.C*d(scsim.Cq.sys.θ,scsim.t,2))) / (2.067833848e-15*scsim.C0.sys.C)
        d(scsim.R2.sys.θ,scsim.t) ~ 3.0385348964360245e15*scsim.R2.sys.R*((3.2910597840193497e-16*d(scsim.C2.sys.θ,scsim.t) - 3.2910597840193497e-16*d(scsim.C0.sys.θ,scsim.t)) / (-scsim.R2.sys.R))
        eqs7[10]
        eqs7[11]
]

sts = [scsim.C1.sys.θ, scsim.C2.sys.θ, scsim.R1.sys.θ, scsim.C0.sys.θ, scsim.R2.sys.θ, scsim.Rq.sys.θ, scsim.Cq.sys.θ]

diff_eq = DifferentialEquation(diff_eqs, sts)

for state in sts
    add_harmonic!(diff_eq, state, scsim.loop1.sys.ω)
end

diff_eq

harmonic_eq = get_harmonic_equations(diff_eq)

#######################################################################################################################
Ntrajectories = 100
tspan = (0.0, 1.0)
tsaves = LinRange(tspan[1]/10.0, tspan[2], 10000)
ps[1] = scsim.loop1.sys.I => 1e-3
ps[8] = scsim.Rq.sys.R => 5.30
ω_vec = ωλ .+ LinRange(-300,300, Ntrajectories)
logger = []
prob = ODAEProblem(model, u0, tspan, ps, saveat = tsaves, maxiters=1e9, force_dtmin=true, dense = false, abstol = 1e-8)
ensemble_prob = EnsembleProblem(prob,prob_func=prob_func, output_func=RMS_IRq)
@time sim2_0p5mA = solve(ensemble_prob,Tsit5(),EnsembleThreads(),trajectories=Ntrajectories)
plot(sim2_0p5mA[2,:])