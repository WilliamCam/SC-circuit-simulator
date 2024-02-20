# In this example we model a DC SQUID driven by a current source

#Load scsim package
include("..//program//SCsim.jl")

#open the circuit netlist file named RF-SQUID.jld2, to create a new circuit netlist we can use scsim.new_netlist()
scsim.open_file("examples//RF-SQUID.jld2")

#cirucit model ODAE system and initial condition vector are created.
model, u0 = scsim.build_circuit()

# we set the values of circuit parameters, for any parameters not specified; default values will be assigned.

I₀ = 1.0e-6
R₀ = 5.0
Φ₀ = scsim.Φ₀

ωc = 2π/Φ₀*I₀*R₀

βc = 0.01
βL = 2*pi/Φ₀ * I₀
βnorm = 2π/Φ₀* I₀ *R₀^2
(8.0/βL)/Φ₀
ps = [
    scsim.loop4.sys.ω => (2*pi*100.0e+6/ωc)
    scsim.loop4.sys.I => 1.25*I₀
    scsim.R1.sys.R => 50.0
    scsim.C1.sys.C => 2.0/βnorm
    scsim.L1.sys.L => 2.0/βL/Φ₀
    scsim.L2.sys.L => 100.0/βL/Φ₀
    scsim.M12.sys.L => 8.0/βL/Φ₀
    scsim.loop1.sys.Φₑ => 0.5
    scsim.β=>βc
    scsim.I₀=>I₀
    scsim.R₀=>R₀
]

tspan = (0.0, 2000.0)
using DifferentialEquations
saveat = LinRange(tspan[2]/10.0, tspan[2], 10000)
prob = ODEProblem(model, u0, tspan, ps, abstol=1e-6, maxiters=1e6)   #Create an ODEProblem to solve for a specified time
sol = solve(prob, ROS3P())

using Plots
plot(sol[scsim.R1.sys.i])

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
