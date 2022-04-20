using DifferentialEquations, Plots, Statistics
include("C://Users//21958742//GitHub//SC-circuit-simulator//program//SCsim.jl")
scsim.open_file("C://Users//21958742//GitHub//SC-circuit-simulator//examples//RLC-nonlinear.jld2")

model, u0, old_model = scsim.build_circuit()

ps = [
    scsim.C1.sys.C => 1.0e-3
    scsim.L1.sys.L => 1.0e-3
    scsim.V1.sys.V => 1.0
    scsim.V1.sys.ω => 1000.0
    scsim.L1.sys.k => 1.0
    scsim.C1.sys.k=>0.0
    scsim.R1.sys.R=>10.0
    scsim.Rin.sys.R=>1.0
]

tspan_ini = (0.0, 1.0)
prob = ODEProblem(model, u0, tspan_ini, ps,  save_everystep = false, progress=true)
sol = solve(prob, Rodas5())

new_u0 = sol[:,end]

N_iniSolves = 50

for i in 1:N_iniSolves
    prob = ODEProblem(model, new_u0, tspan_ini, ps,  save_everystep = false)
    sol = solve(prob, Rodas5())
    new_u0 = sol[:,end]
end

tspan = (0.0, 5.0)

tsaves = LinRange(3.0, tspan[2], 10000)
dt = tsaves[2] - tsaves[1]


prob = ODEProblem(model, new_u0, tspan, ps, saveat = tsaves, maxiters=1e8)
@time sol = solve(prob, ROS3P();)

v =  1/dt * scsim.Φ₀/(2*pi)*diff(sol[scsim.R1.sys.θ])
i = sol[scsim.Rin.sys.i][2:end]
plot(v)
fs = 2/dt
using FFTW
v_fft = 2*rfft(v)/length(v)
fn = LinRange(0.0, fs/2, length(v_fft))
plot(fn, abs.(v_fft))
sum(abs.(v_fft))

function RMS(x)
    return sqrt(mean((x .- mean(x)).^2))
end
RMS(v)

Ntrajectories = 500

ωλ = 1/sqrt(prob.p[4]*prob.p[6])
ω_vec = ωλ .+ LinRange(-0.6*ωλ, 0.6*ωλ, Ntrajectories)

function prob_func(prob, i ,repeat)
    prob.p[13] = ω_vec[i]
    prob
end


function output_func(sol,i)
    (RMS(sol[scsim.Rin.sys.i]),false)
end

ensemble_prob = EnsembleProblem(prob,prob_func=prob_func, output_func=output_func)
@time sim = solve(ensemble_prob,Rodas5(),EnsembleThreads(),trajectories=Ntrajectories)

plot(ω_vec/(2*pi), sim.u)

ps = [
    scsim.C1.sys.C => 1.0e-3
    scsim.L1.sys.L => 1.0e-3
    scsim.V1.sys.V => 50.0
    scsim.V1.sys.ω => 1000.0
    scsim.L1.sys.k => 1.0
    scsim.C1.sys.k=>0.0
    scsim.R1.sys.R=>10.0
    scsim.Rin.sys.R=>1.0
]

prob = ODEProblem(model, new_u0, tspan, ps, saveat = tsaves, maxiters=1e8)
ensemble_prob = EnsembleProblem(prob,prob_func=prob_func, output_func=output_func)
@time sim2 = solve(ensemble_prob,Rodas5(),EnsembleThreads(),trajectories=Ntrajectories)
ps = [
    scsim.C1.sys.C => 1.0e-3
    scsim.L1.sys.L => 1.0e-3
    scsim.V1.sys.V => 1e+3
    scsim.V1.sys.ω => 1000.0
    scsim.L1.sys.k => 1.0
    scsim.C1.sys.k=>0.0
    scsim.R1.sys.R=>10.0
    scsim.Rin.sys.R=>1.0
]

prob = ODEProblem(model, new_u0, tspan, ps, saveat = tsaves, maxiters=1e8)
ensemble_prob = EnsembleProblem(prob,prob_func=prob_func, output_func=output_func)
@time sim3 = solve(ensemble_prob,ROS3P(),EnsembleSerial(),trajectories=Ntrajectories)

ps = [
    scsim.C1.sys.C => 1.0e-3
    scsim.L1.sys.L => 1.0e-3
    scsim.V1.sys.V => 2e+3
    scsim.L1.sys.k => 1.0
    scsim.C1.sys.k=>0.0
    scsim.R1.sys.R=>10.0
    scsim.Rin.sys.R=>1.0
    scsim.Rin.sys.k=>0.0
    scsim.R1.sys.k=>0.0
    scsim.loop1.sys.Φₑ=>0.0
    scsim.loop2.sys.Φₑ=>0.0
    scsim.loop3.sys.Φₑ=>0.0
]

prob = ODEProblem(model, new_u0, tspan, ps, saveat = tsaves, maxiters=1e8)
ensemble_prob = EnsembleProblem(prob,prob_func=prob_func, output_func=output_func)
@time sim4 = solve(ensemble_prob,Rodas5(),EnsembleThreads(),trajectories=Ntrajectories)

using HarmonicBalance, ModelingToolkit
eqs = equations(model)
sts = states(model)

#eqs = substitute(eqs, Dict(ps))

eqs1 = [scsim.D2(scsim.C1.sys.θ) ~ (6.283185307179586*scsim.C1.sys.i) / 2.0678338480000003e-18]
eqs1 = substitute(eqs1, Dict(scsim.C1.sys.i => scsim.Rin.sys.i - scsim.loop2.sys.iₘ))
append!(eqs1,eqs[3:end])
eqs1 = substitute(eqs1, Dict(sts[1]=>scsim.D(scsim.C1.sys.θ)))

diff_eq = DifferentialEquation(eqs1, [scsim.C1.sys.θ, scsim.R1.sys.θ, scsim.Rin.sys.θ, 
scsim.V1.sys.θ, scsim.R1.sys.i, scsim.Rin.sys.i, scsim.loop2.sys.iₘ])

add_harmonic!(diff_eq, scsim.C1.sys.θ, scsim.V1.sys.ω)
add_harmonic!(diff_eq, scsim.R1.sys.θ, scsim.V1.sys.ω)
add_harmonic!(diff_eq, scsim.Rin.sys.θ, scsim.V1.sys.ω)
add_harmonic!(diff_eq, scsim.V1.sys.θ, scsim.V1.sys.ω)
add_harmonic!(diff_eq, scsim.R1.sys.i, scsim.V1.sys.ω)
add_harmonic!(diff_eq, scsim.Rin.sys.i, scsim.V1.sys.ω)
add_harmonic!(diff_eq, scsim.loop2.sys.iₘ, scsim.V1.sys.ω)

harmonic_eq = get_harmonic_equations(diff_eq)

varied = scsim.V1.sys.ω => LinRange(0.9e+3, 1.2e+3, 100)
fixed = (
    scsim.C1.sys.C => 1.0e-3,
    scsim.L1.sys.L => 1.0e-3,
    scsim.V1.sys.V => 2e+3,
    scsim.L1.sys.k => 1.0,
    scsim.C1.sys.k=>1.0,
    scsim.R1.sys.R=>10.0,
    scsim.Rin.sys.R=>1.0,
    scsim.Rin.sys.k=>1.0,
    scsim.R1.sys.k=>1.0,
    scsim.loop1.sys.Φₑ=>1.0,
    scsim.loop2.sys.Φₑ=>1.0,
    scsim.loop3.sys.Φₑ=>1.0
)

solutions = get_steady_states(harmonic_eq, varied, fixed)