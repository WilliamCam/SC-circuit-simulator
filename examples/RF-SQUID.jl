# In this example we model a DC SQUID driven by a current source

#Load scsim package
include("C://Users//21958742//GitHub//SC-circuit-simulator//program//SCsim.jl")

#open the circuit netlist file named RF-SQUID.jld2, to create a new circuit netlist we can use scsim.new_netlist()
scsim.open_file("C://Users//21958742//GitHub//SC-circuit-simulator//examples//RF-SQUID.jld2")

#cirucit model ODAE system and initial condition vector are created.
model, u0 = scsim.build_circuit()

# we set the values of circuit parameters, for any parameters not specified; default values will be assigned.

I₀ = 1.0e-6
R₀ = 5.0
Φ₀ = scsim.Φ₀

βc  = 2*pi/Φ₀ * I₀ * R₀^2
βL = 2*pi/Φ₀ * I₀

ps = [
    scsim.loop4.sys.ω => 2*pi*100.0e+6
    scsim.loop4.sys.I => 1.25*I₀
    scsim.J1.sys.I0 => I₀
    scsim.J1.sys.R => R₀
    scsim.J1.sys.C => 0.01/βc
    scsim.J1.sys.L => 1.0/βc
    scsim.R1.sys.R => 50.0
    scsim.C1.sys.C => 2.0/βc
    scsim.L1.sys.L => 2.0/βL
    scsim.L2.sys.L => 100.0/βL
    scsim.M12.sys.L => 8.0/βL
    scsim.loop1.sys.Φₑ => 0.5*Φ₀
]

tspan = (0.0, 1e-6)
saveat = LinRange(tspan[2]/10.0, tspan[2], 10000)

sol = scsim.tsolve(model, u0, tspan, ps, saveat = saveat)

scsim.tplot(sol, scsim.R1, units="volts")


Φspan = (0.0, 2.0*Φ₀)

ensemble_sol = scsim.ensemble_parameter_sweep(
    model, u0, tspan, Φspan, ps, scsim.loop1.sys.Φₑ, scsim.R1, saveat = saveat
)
using Plots
plot(ensemble_sol.u)

Ispan = (0.0, 2*I₀)
ensemble_sol = scsim.ensemble_parameter_sweep(
    model, u0, tspan, Ispan, ps, scsim.loop4.sys.I, scsim.R1, saveat = saveat
)
plot(ensemble_sol.u)
