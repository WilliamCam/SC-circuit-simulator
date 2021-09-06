module Simulate

using DifferentialEquations
using Statistics
using Plots
include("constants.jl")
include("EquationsOfMotion.jl")
include("Sweeps\\FrequencySweep.jl")

"""
    RMS(x)

Computes the root mean squared of a vector x
"""
function RMS(x)
    sqrt(mean(x.^2))
end
"""
    output_func(sol,i)
output function for ensemble problem.
Gives the RMS Voltage.
"""
function output_func(sol, i)
    (RMS(@. phi0/(2.0*pi)*sol[1,:]),false)
end

function Frequency_response(Ib, Iin, NΦ₀, Nperiods, N_f_values, df)
    tspan = (0.0,2*pi*Nperiods/ωλ)
    tsaves = LinRange(0.0,tspan[end],1000);

    df_vec = LinRange(-df*ωλ,df*ωλ,N_f_values)

    p_consts = df_sweep_consts(Ib, Iin, NΦ₀*phi0)
    p_vars = df_sweep_vars(ωλ+df_vec[1],junc_i)

    p = df_sweep_params(p_consts,p_vars);

    function prob_func(prob,i,repeat)
        prob.p.v.ω = ωλ + df_vec[i]
        prob
    end

    #solve for inital conditions
    prob_RLC_para_ini = ODEProblem(RLC_para_fswp!, u0, (0.0, 100/ωλ), p, abstol = 1e-3);
    ui = solve(prob_RLC_para_ini, save_everystep=false).u[end];

    prob_RLC_para = ODEProblem(RLC_para_fswp!, ui, tspan, p, abstol = 1e-3, saveat = tsaves, save_idxs=[10]);

    ensemble_prob = EnsembleProblem(prob_RLC_para, prob_func=prob_func, output_func = output_func);

    sim = solve(ensemble_prob,Tsit5(),EnsembleThreads(),trajectories=N_f_values, maxiters=1e9, progress=true)

    graph = plot((ωλ .+ df_vec)./ωλ,sim[:],
        title = "RLC-SQUID Parallel Config \n Frequency Response", xlabel = "Detuning (δ)", ylabel = "Vout (Vrms)")

    display(graph)
    return sim
end
export Frequency_response
end
