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
include("Sweeps\\flux.jl")

"""
    single_sim(Ib::float64, Iin::float64, NΦe::float64; ωin = ωλ, tsim = 5*2.0*pi/ωλ)

Performs a single transient simulation of the RLC-SQUID system with bias current Ib, input current Ib*sin(ωin*t)
flux bias NΦe*Φ₀ to t = tsim.
"""
function single_sim(Ib::Float64, Iin::Float64, NΦe::Float64; ωin = ωλ, tsim = 3*2.0*pi/ωλ)
    junc_i = zeros(5) #junction current cache initialisation
    u0 = zeros(11) #junction phase cache initialisation

    tspan = (0.0, tsim) #simulation timespan
    ini_tspan = (0.0, 1e-6) #inital value simulation timespan
    tsaves = LinRange(2.0*pi/ωλ,tspan[end],10000) #points to save at

    p_consts = flux_sweep_consts(Ib, Iin, ωin) #input constant parameters
    p_vars = flux_sweep_vars(NΦe*phi0,junc_i) #input variable parameters
    p = flux_sweep_params(p_consts,p_vars); #input parameter struct

    ini_prob = ODEProblem(RLC_para_v2_fluxswp!, u0, ini_tspan, p, abstol = 1e-5, reltol=1e-5) #initial problem
    ui = solve(ini_prob, Tsit5(), maxiters=1e9, progress=true, save_everystep=false)[:,end] #solve for initial conditions

    prob_RLC_para = ODEProblem(RLC_para_v2_fluxswp!, ui, tspan, p, abstol = 1e-5, reltol=1e-5, saveat = tsaves); #problem

    sim = solve(prob_RLC_para, Tsit5(), maxiters=1e9, progress=true) #solve problem

    V_input = phi0/(2*pi)*sim[6,:]
    I_input = @. Iin*sin(ωin*sim.t) - V_input*Ginp

    V_out = phi0/(2*pi)*sim[8,:]
    I_out = phi0/(2*pi)*C2*sim[11,:]

    graph = plot(sim.t, [V_out,V_input],
        title = "RLC-SQUID Parallel Config",
        xlabel = "Time (s)",
        ylabel = "Voltage (V)"
    )
    display(graph)


    return sim
end #single_sim


end #end module
