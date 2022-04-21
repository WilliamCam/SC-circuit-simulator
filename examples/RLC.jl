include("C://Users//21958742//GitHub//SC-circuit-simulator//program//SCsim.jl")
scsim.open_file("C://Users//21958742//GitHub//SC-circuit-simulator//examples//RLC.jld2")
model, u0 = scsim.build_circuit()

ps = []
scsim.parameter_set(ps, scsim.C1.sys.C, 1.0)
scsim.parameter_set(ps, scsim.V1.sys.V, 1.0)
scsim.parameter_set(ps, scsim.V1.sys.ω, 1.0)

u_initial = scsim.solve_ini(model, u0, 1.0, ps)

prob = ODEProblem(model, u0, (0.0,1.0), ps)

solve(prob)

tspan = (0.0, 20.0)

sol = scsim.tsolve(model, u_initial, tspan, ps)

scsim.tplot(sol, scsim.R1)

u0