module Simulate
#simulation module for RLC-SQUID Parallel Configuration with low pass filter

using DifferentialEquations
using Statistics
using Plots
using FileIO

include("constants.jl")
include("EquationsOfMotion.jl")
include("Sweeps\\FrequencySweep.jl")

"""
    RMS(x)

Computes the root mean squared of a vector x
"""
function RMS(x)
    sqrt(mean(x.^2))
end #end RMS

"""
    output_func(sol,i)

output function for ensemble problem.
Gives the RMS Voltage.
"""
function output_func(sol, i)
    (RMS(@. phi0/(2.0*pi)*sol[1,:]),false)
end #end output_func

"""
    frequency_response(Ib, Iin, NΦ₀, df, N_f_values = 500, N_periods=50

Simulates the frequency response of the system from -df*ωλ to df*ωλ by
performing N_f_values transient simulations 2*pi*N_periods/ωλ seconds long
"""
function frequency_response(Ib, Iin, NΦ₀, df; N_fvalues = 500, N_periods=100, export_result = false)
    tspan = (0.0,2*pi*N_periods/ωλ) #simulation timespan
    tsaves = LinRange(2*pi/ωλ,tspan[end],1000); #time points to save at

    df_vec = LinRange(-df*ωλ,df*ωλ,N_fvalues) #frequency points

    p_consts = df_sweep_consts(Ib, Iin, NΦ₀*phi0) #input constant parameters
    p_vars = df_sweep_vars(ωλ+df_vec[1],junc_i) #input variable parameters

    p = df_sweep_params(p_consts,p_vars); #input parameters master struct

    function prob_func(prob,i,repeat) #problem funtion modifies input parameters
        prob.p.v.ω = ωλ + df_vec[i]
        prob
    end #end prob_func

    prob_RLC_para = ODEProblem(RLC_para_fswp!, u0, tspan, p, abstol = 1e-3, saveat = tsaves, save_idxs=[10]); #problem

    ensemble_prob = EnsembleProblem(prob_RLC_para, prob_func=prob_func, output_func = output_func); #ensemble

    sim = solve(ensemble_prob,Vern6(),EnsembleThreads(),trajectories=N_fvalues, maxiters=1e9, progress=true) #solve ensemble

    #plot result
    graph = plot((ωλ .+ df_vec)./ωλ,sim[:],
        title = "RLC-SQUID Parallel Config \n Frequency Response",
        xlabel = "Detuning (δ)",
        ylabel = "Vout (Vrms)"
    )
    display(graph)

    if export_result == true #export simulation result to .jld2
        label_string = string(p_consts.Iin*1.0e9) * "nV-" * string(NΦ₀) * "Φ₀-" * string(df) * "δ"
        with jldopen("simlogs.jld2", "a+") do file
            file[label_string] = sim
    end #data export

    return sim


end #end frequency_response

export frequency_response

end #end module
