#In this example we model an parallel RLC circuit driven by a voltage source

#Load scsim package
include("C://Users//21958742//GitHub//SC-circuit-simulator//program//SCsim.jl")

#open the circuit netlist file named RLC.jld2, to create a new circuit netlist we can use scsim.new_netlist()
scsim.open_file("C://Users//21958742//GitHub//SC-circuit-simulator//examples//RLC.jld2")

#cirucit model ODE system and initial condition vector are created.
model, u0 = scsim.build_circuit()

# we set the values of circuit parameters, for any parameters not specified; default values will be assigned.
ps = [
    scsim.V1.sys.ω => 1.0
    scsim.V1.sys.V => 1.0
    scsim.C1.sys.C => 1.0
]

#find resaonable initial conditions by settling transients
u_initial = scsim.solve_ini(model, u0, 1.0, ps)

#specify transient window for solver
tspan = (0.0, 30.0)

#transient circuit analysis
sol = scsim.tsolve(model, u_initial, tspan, ps)
scsim.tplot(sol, scsim.R1, units = "Volts")

# we can pass any arguments known to the problem interface of DifferentialEquations.jl
using DifferentialEquations

tsaves = LinRange(10.0, 30.0, 500)
sol = scsim.tsolve(model, u_initial, tspan, ps, saveat = tsaves, alg = DifferentialEquations.ROS3P())
scsim.tplot(sol, scsim.R1, units = "Volts")
