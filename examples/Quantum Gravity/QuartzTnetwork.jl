using DifferentialEquations, Plots, Statistics
include("C://Users//21958742//GitHub//SC-circuit-simulator//program//SCsim.jl")
scsim.open_file("C://Users//21958742//GitHub//SC-circuit-simulator//examples//QuartzTnetwork.jld2")

model, u0, old_model = scsim.build_circuit()

ps = [
    scsim.V1.sys.V => 1.0e-6
    scsim.V1.sys.ω => 5.0e+7
    scsim.R1.sys.R=> 50.0
    scsim.C1.sys.C => 10.0e-9
    scsim.Cq.sys.C => 1.0e-16
    scsim.Lq.sys.L => 4.0
    scsim.Lq.sys.k => 1.0
    scsim.Rq.sys.R => 100.0
    scsim.C0.sys.C=> 1.0e-11
    scsim.C2.sys.C => 10.0e-9
    scsim.R2.sys.R=>1.0e+3
]
tspan_ini = (0.0, 1.0e-6)
prob = ODEProblem(model, u0, tspan_ini, ps,  save_everystep = false, progress=true)
sol = solve(prob, ROS3P())

new_u0 = sol[:,end]


tspan = (0.0, 1.0e-3)
tsaves = LinRange(tspan[2]/10.0, tspan[2], 10000)
dt = tsaves[2] - tsaves[1]


prob = ODEProblem(model, u0, tspan, ps, saveat = tsaves, maxiters=1e9, abstol = 1e-12)
@time sol = solve(prob, ROS3P())
v =  1/dt * scsim.Φ₀/(2*pi)*diff(sol[scsim.R2.sys.θ])
i = sol[scsim.R2.sys.i][2:end]
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

Ntrajectories = 50

ωλ = 1/sqrt(prob.p[12]*prob.p[2])
ω_vec = ωλ .+ LinRange(-0.01*ωλ, 0.01*ωλ, Ntrajectories)

function prob_func(prob, i ,repeat)
    prob.p[7] = ω_vec[i]
    prob
end


function output_func(sol,i)
    (RMS(sol[scsim.R2.sys.i]),false)
end

ensemble_prob = EnsembleProblem(prob,prob_func=prob_func, output_func=output_func)
@time sim = solve(ensemble_prob,ROS3P(),EnsembleThreads(),trajectories=Ntrajectories)

plot(ω_vec/(2*pi), sim.u)