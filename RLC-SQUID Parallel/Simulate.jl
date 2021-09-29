
module Simulate
#simulation module for RLC-SQUID Parallel Configuration with low pass filter

using DifferentialEquations
using Statistics
using Plots
using JLD2
using FileIO

include("constants.jl")
include("EquationsOfMotion.jl")
include("Sweeps\\FrequencySweep.jl")
include("Sweeps\\FluxSweep.jl")

"""
    RMS(x)

Computes the root mean squared of a vector x after removing any DC offset
"""
function RMS(x)
    sqrt(mean((x.-mean(x)).^2))
end #end RMS

"""
    output_func(sol,i)

output function for ensemble problem.
Gives the RMS Voltage.
"""
function output_func_Vrms(sol, i)
    (RMS(@. phi0/(2.0*pi)*sol[1,:]),false)
end #end output_func

"""
    frequency_response(Ib, Iin, NΦ₀, df; N_fvalues = 500, N_periods=100, export_result = false)

Simulates the frequency response of the system from -df*ωλ to df*ωλ by
performing N_f_values parallel transient simulations 2*pi*N_periods/ωλ seconds long
and returning the average RMS voltage at each frequency
"""
function freq_response(Ib, Iin, NΦ₀, df; N_fvalues = 500, N_periods=100, export_result = false)
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

    ensemble_prob = EnsembleProblem(prob_RLC_para, prob_func=prob_func, output_func = output_func_Vrms); #ensemble

    sim = solve(ensemble_prob,Vern6(),EnsembleThreads(),trajectories=N_fvalues, maxiters=1e9, progress=true) #solve ensemble

    #plot result
    graph = plot((ωλ .+ df_vec)./ωλ,sim[:],
        title = "RLC-SQUID Parallel Config \n Frequency Response",
        xlabel = "Detuning (δ)",
        ylabel = "Vout (Vrms)"
    )
    display(graph)

    if export_result == true #export simulation result to .jld2
        label_string = "FreqResp-" * string(p_consts.Iin*1.0e9) * "nV-" * string(NΦ₀) * "Φ₀-" * string(df) * "δ"
        jldopen("simlogs.jld2", "a+") do file
            file[label_string] = sim
        end # file open

    end #data export

    return sim
end #end frequency_response
"""
    function S21_single_sim(Ib::Float64, Iin::Float64, NΦ₀::Float64; N_periods = 100, ωin = ωλ)

Performs a single transient simulation using FrequencySweep.jl in order to determine SQUID amplifier
gain.
"""
function S21_single_sim(Ib::Float64, Iin::Float64, NΦ₀::Float64; N_periods = 100, ωin = ωλ)
    tspan = (0.0,2*pi*N_periods/ωλ) #simulation timespan
    tsaves = LinRange(2*pi/ωλ,tspan[end],10000); #time points to save at

    p_consts = df_sweep_consts(Ib, Iin, NΦ₀*phi0) #input constant parameters
    p_vars = df_sweep_vars(ωin,junc_i) #input variable parameters
    p = df_sweep_params(p_consts,p_vars); #input parameters master struct

    prob_RLC_para = ODEProblem(RLC_para_fswp_S21!, u0_S21, tspan, p, abstol = 1e-3, saveat = tsaves, save_idxs=[10,11]); #problem

    # saved_values = SavedValues(Float64, Array{Float64})
    # cb = SavingCallback((u, t, integrator)->integrator(t,Val{0}), saved_values, saveat = tsaves)

    sim = solve(prob_RLC_para,Tsit5(), maxiters=1e9, progress=true) #solve problem

    Vin = Rλ * Iin .* sin.((ωλ).*sim.t) #Total available power on resonance
    Vout = phi0/(2.0*pi).*sim[1,:]

    Iout = phi0/(2.0*pi)*C2.*sim[2,:]

    Vgain = RMS(Vout)/RMS(Vin)
    Igain = RMS(Iout)/(0.707*Iin)
    Pgain = Vgain*Igain



    graph = plot(sim.t, [Vout.-mean(Vout)],
        title = "RLC-SQUID Parallel Config \n S21 parameter sim",
        xlabel = "Time (s)",
        ylabel = "Voltage (V)"
    )
    display(graph)

    return  RMS(Vout), Igain, Pgain
end #end S21_single_sim
"""
    flux_response(Ib::Float64, Iin::Float64; N_periods = 100, ωin = ωλ, Npoints = 200)

Determines the flux-voltage response VΦ of the RLC-SQUID system by perfroming NPoints
parallel transient simulations with varying flux bias from 0 to 2Φ₀, solving for
the mean rms output voltage of the system. Result is plotted and saved to simlogs.jld2
if export_result = true.
"""
function flux_response(Ib::Float64, Iin::Float64; N_periods = 20, ωin = ωλ, Npoints = 200, export_result=true)
    tspan = (0.0,2*pi*N_periods/ωλ) #simulation timespan
    tsaves = LinRange(2*pi/ωλ,tspan[end],1000); #time points to save at

    flux_vec = LinRange(0.0,2.0*phi0,Npoints) #frequency points

    p_consts = flux_sweep_consts(Ib, Iin, ωin) #input constant parameters
    p_vars = flux_sweep_vars(flux_vec[1], junc_i) #input variable parameters
    p = flux_sweep_params(p_consts,p_vars); #input parameters master struct

    function prob_func(prob,i,repeat) #problem funtion modifies input parameters
        prob.p.v.Φe = flux_vec[i]
        prob
    end #end prob_func

    function output_func(sol, i)
        (mean(phi0/(2.0*pi)*1.0e+6*sol),false)
    end# output_func

    prob = ODEProblem(RLC_para_fluxswp!, u0, tspan, p, abstol = 1e-5, saveat = tsaves, save_idxs=[10]); #problem

    ensemble_prob = EnsembleProblem(prob, prob_func=prob_func, output_func = output_func); #ensemble

    sim = solve(ensemble_prob,Tsit5(),EnsembleThreads(),trajectories=Npoints, maxiters=1e9, progress=true) #solve ensemble

    #plot result
    graph = plot(flux_vec/phi0,sim[:],
        title = "RLC-SQUID Parallel Config \n Flux-Voltage Response",
        xlabel = "Flux (Φ₀)",
        ylabel = "Averaged Vout (μVrms)"
    )
    display(graph)

    if export_result == true #export simulation result to .jld2
        label_string = "FluxResp-" * string(p_consts.Iin*1.0e9) * "nV-" * string(ωin/ωλ) * "ωλ"
        jldopen("simlogs.jld2", "a+") do file
            file[label_string] = sim
        end # file open

    end #data export

    return sim
end #end flux_response


####end module
export frequency_response
end #end module
