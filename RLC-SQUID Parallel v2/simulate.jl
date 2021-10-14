module simulate
#simulation module for RLC-SQUID Parallel configuration with low pass filter v2
#second version includes parasitic capacitacne C₀ due to BAW electrodes, as well as input source impedance Rinp

using DifferentialEquations
using Statistics
using Plots
using JLD2
using FileIO

include("constants.jl")
include("system.jl")
include("models.jl")

"""
    initial_state_sim(Ib::Float64, Iin::Float64, NΦe::Float64; ωin = ωλ, tsim = 1.0e-6)

Performs a transient simulation in order to find a final state vector u0 to be used as an initial
state for further simulations, this helps settle initial simulation transients.
returns the final state vector after time tsim; ui
"""
function initial_state_sim(Ib::Float64, Iin::Float64, NΦe::Float64; ωin = ωλ, tsim = 1.0e-6)
    u0 = zeros(11) #state vector cache initialisation
    junc_i = zeros(5) #junction current cache initialisation

    tspan = (0.0, tsim) #simulation timespan

    p_consts = flux_sweep_consts(Ib, Iin, ωin) #input constant parameters
    p_vars = flux_sweep_vars(NΦe*phi0,junc_i) #input variable parameters
    p = flux_sweep_params(p_consts,p_vars); #input parameter struct

    prob = ODEProblem(RLC_para_v2_fluxswp!, u0, tspan, p, abstol = 1e-4, reltol=1e-6) #initial problem
    ui = solve(prob, Tsit5(), maxiters=1e9, progress=true, save_everystep=false)[:,end] #solve for initial conditions

    return ui #final state vector
end #end initial_state_sim
"""
single_sim(u0::Vector{Float64}, Ib::Float64, Iin::Float64, NΦe::Float64;
     ωin = ωλ, tsim = 3*2.0*pi/ωλ)

Performs a single transient simulation of the RLC-SQUID system with bias current
Ib, input current Ib*sin(ωin*t) flux bias NΦe*Φ₀ to t = tsim.
returns solution object sol with time vector sol.t and solution sol.u
"""
function single_sim(u0::Vector{Float64}, Ib::Float64, Iin::Float64, NΦe::Float64;
     ωin = ωλ, tsim = 10*2.0*pi/ωλ, tstart = 0.0
     )
     junc_i = zeros(5) #junction current cache initialisation

    tspan = (0.0, tsim) #simulation timespan
    tsaves = LinRange(tstart,tspan[end],50000) #points to save at

    p_consts = flux_sweep_consts(Ib, Iin, ωin) #input constant parameters
    p_vars = flux_sweep_vars(NΦe*phi0,junc_i) #input variable parameters
    p = flux_sweep_params(p_consts,p_vars); #input parameter struct


    prob_RLC_para = ODEProblem(RLC_para_v2_fluxswp!, u0, tspan, p,
        abstol = 1e-4, reltol = 1e-6, saveat = tsaves
    ); #problem

    sim = solve(prob_RLC_para, Tsit5(), maxiters=1e9, progress=true) #solve problem

    V_input = phi0/(2*pi)*sim[6,:]
    I_input = @. Iin*sin(ωin*sim.t) - V_input * Ginp

    V_out = phi0/(2*pi)*sim[10,:]
    I_out = phi0/(2*pi)*C2*sim[11,:]

    graph = plot(sim.t/τ, [V_input*1e+6, (V_out.-mean(V_out))*1e+6],
        title = "RLC-SQUID Parallel Config",
        xlabel = "Time (τ)",
        ylabel = "Voltage (V)"
    )
    display(graph)

    return sim, V_input, V_out, I_input, I_out

end #single_sim

"""
    flux_voltage(Ib::Float64, Iin::Float64; ωin = ωλ, tsim = 10*2*pi/ωλ, Npts = 200, export_result = true)

Performs Npts transient simulations of length tsim for varying flux bias input. Returns
the external flux - voltage relationship of the system. result can be saved to simlogs.jl2 if
export_result = true.

"""
function flux_voltage(Ib::Float64, Iin::Float64; ωin = ωλ, tsim = 10*2*pi/ωλ, Npts = 200, export_result = true)
    u0 = zeros(11) #state vector cache initialisation
    junc_i = zeros(5) #junction current cache initialisation

    flux_vec = phi0*LinRange(0.0, 2.0, Npts) #flux bias values

    tspan = (0.0, tsim) #simulation timespan
    ini_tspan = (0.0, 1e-6) #inital value simulation timespan
    tsaves = LinRange(0.0,tspan[end],10000) #points to save at

    p_consts = flux_sweep_consts(Ib, Iin, ωin) #input constant parameters
    p_vars = flux_sweep_vars(flux_vec[1],junc_i) #input variable parameters
    p = flux_sweep_params(p_consts,p_vars); #input parameter struct



    function prob_func(prob,i,repeat) #problem funtion modifies input parameters
        prob.p.v.Φe = flux_vec[i]
        prob
    end #end prob_func

    function output_func(sol, i)   #outputs the mean voltage across C₂
        (mean(phi0/(2.0*pi)*sol),false)
    end# output_func

    prob = ODEProblem(RLC_para_v2_fluxswp!, u0, tspan, p, abstol = 1e-5,saveat = tsaves, save_idxs=[10]); #problem

    ensemble_prob = EnsembleProblem(prob, prob_func=prob_func, output_func = output_func); #ensemble

    sim = solve(ensemble_prob,Tsit5(),EnsembleThreads(),trajectories=Npts, maxiters=1e9, progress=true) #solve ensemble

    #plot result
    graph = plot(flux_vec/phi0,1.0e+6.*sim[:],
        title = "RLC-SQUID Parallel Config v2 \n Flux-Voltage Response",
        xlabel = "Flux (Φ₀)",
        ylabel = "Averaged Vout (μVrms)"
    )
    display(graph)

    if export_result == true #export simulation result to .jld2
        label_string = "VΦ-" * string(p_consts.Iin*1.0e9) * "nV-" * string(ωin/ωλ) * "ωλ"
        jldopen("simlogs.jld2", "a+") do file
            file["RLCPara-v2/" * label_string] = sim
        end # file open

    end #data export
    return sim
end #flux_voltage

function S21_sim(Ib::Float64, Iin::Float64, Φe::Float64;
        δ = 0.5, tsim = 10*2*pi/ωλ, Npts = 200, export_result = true,
        tstart = 3*2*pi/ωλ, S21 = true
    )

    u0 = zeros(11) #state vector cache initialisation
    junc_i = zeros(5) #junction current cache initialisation

    f_vec = ωλ .+ ωλ.*LinRange(-δ, δ, Npts) #input frequency vector

    tspan = (0.0, tsim) #simulation timespan
    tsaves = LinRange(tstart,tspan[end],10000) #points to save at

    p_consts = freq_sweep_consts(Ib, Iin, Φe) #input constant parameters
    p_vars = freq_sweep_vars(f_vec[1],junc_i) #input variable parameters
    p = freq_sweep_params(p_consts,p_vars); #input parameter struct

    function mean_RMS(sol)
        sqrt(mean((sol .- mean((sol))).^2))
    end # end mean_RMS

    function prob_func(prob,i,repeat) #problem funtion modifies input parameters
        prob.p.v.ω = f_vec[i]
        prob
    end #end prob_func



    function S21_func(sol, i)   #outputs the mean RMS voltage across C₂
        ((mean_RMS(phi0/(2*pi)*sol[2,:])*mean_RMS(C2*phi0/(2*pi)*sol[3,:]))/(mean_RMS(phi0/(2*pi)*sol[1,:])*mean_RMS(Iin* sin.(f_vec[i]*sol.t) .- phi0/(2*pi)*sol[1,:] .* Ginp)),false)
    end# end S21

    function VoutRMS_func(sol, i)   #outputs the mean RMS voltage across C₂
        (mean_RMS(phi0/(2*pi)*sol[2,:]),false)
    end# end S21

    if S21 == true
        output_func = S21_func
    else
        output_func = VoutRMS_func
    end #end conditional

    prob = ODEProblem(RLC_para_v2_fswp!, u0, tspan, p, abstol = 1e-5,saveat = tsaves, save_idxs=[6,10,11]); #problem

    ensemble_prob = EnsembleProblem(prob, prob_func=prob_func, output_func = output_func); #ensemble

    sim = solve(ensemble_prob,Tsit5(),EnsembleThreads(),trajectories=Npts, maxiters=1e9, progress=true) #solve ensemble

    #plot result
    graph = plot(f_vec/(2*pi),sim[:],
        title = "RLC-SQUID Parallel Config v2 \n Amplitude Response",
        xlabel = "Frequency (Hz)",
        ylabel = "Averaged Vout (μVrms)"
    )
    display(graph)

    if export_result == true #export simulation result to .jld2
        label_string = "S21-" * string(p_consts.Iin*1.0e9) * "nV"
        jldopen("simlogs.jld2", "a+") do file
            file["RLCPara-v2/" * label_string] = sim
        end # file open

    end #data export
    return sim
end #end S21_sim



end #end module
