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
    scsim.J1.sys.L => 1.0/βL
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

ωspan = (0.1*2*pi*100.0e+6, 2*2*pi*100.0e+6)
ensemble_sol = scsim.ensemble_parameter_sweep(
    model, u0, tspan, ωspan, ps, scsim.loop4.sys.ω, scsim.R1, saveat = saveat
)
plot(ensemble_sol.u)

#### Using Harmonoic Balance
using HarmonicBalance


using ModelingToolkit

obs = observed(model)
eqs = equations(model)

obs[6]


function is_term(eqn, target_term)
    vars = get_variables(eqn)
    ret = false
    for term in vars
        if isequal(term, target_term)
            ret = true
            break
        else
            ret = false
        end
    end
    return ret
end



loop2_eqs = [eqn for eqn in obs if is_term(eqn, scsim.loop2.sys.iₘ)]

eqs2 = substitute(eqs, Dict(scsim.loop2.sys.iₘ => Symbolics.solve_for(loop2_eqs[1], scsim.loop2.sys.iₘ)))

J1_eqs = [eqn for eqn in obs if is_term(eqn, scsim.J1.sys.i)]

obs[6]

eqs3 = substitute(obs[6], Dict(scsim.J1.sys.i => Symbolics.solve_for(J1_eqs[3], scsim.J1.sys.i)))

eqs4 = substitute(eqs3, Dict(scsim.loop2.sys.iₘ => Symbolics.solve_for(loop2_eqs[2], scsim.loop2.sys.iₘ)))

C1_eqs = [eqn for eqn in obs if is_term(eqn, scsim.C1.sys.i)]

eqs5 = substitute(eqs4, Dict(scsim.C1.sys.i => Symbolics.solve_for(C1_eqs[2], scsim.C1.sys.i)))

loop3_eqs = [eqn for eqn in obs if is_term(eqn, scsim.loop3.sys.iₘ)]

eqs6 = substitute(eqs5, Dict(scsim.loop3.sys.iₘ => Symbolics.solve_for(loop3_eqs[1], scsim.loop3.sys.iₘ)))

R1_eqs = [eqn for eqn in obs if is_term(eqn, scsim.R1.sys.i)]

eqs7 = substitute(eqs6, Dict(scsim.R1.sys.i => Symbolics.solve_for(R1_eqs[1], scsim.R1.sys.i)))

dRθ = get_variables(obs[4])[1]

eqs8 = substitute(eqs7, Dict(dRθ => Symbolics.solve_for(obs[4], dRθ)))

loop4_eqs =  [eqn for eqn in obs if is_term(eqn, scsim.loop4.sys.iₘ)]

eqs9 = substitute(eqs8, Dict(scsim.loop4.sys.iₘ => Symbolics.solve_for(loop4_eqs[1], scsim.loop4.sys.iₘ)))

eqs10 = Symbolics.simplify(eqs9)

eqs11 = obs[12]

eqs12 = substitute(eqs11, Dict(scsim.C1.sys.i => Symbolics.solve_for(C1_eqs[1], scsim.C1.sys.i)))

eqs13= substitute(eqs12, Dict(scsim.loop2.sys.iₘ => Symbolics.solve_for(eqs[5], scsim.loop2.sys.iₘ)))

eqs14 = substitute(eqs13, Dict(scsim.J1.sys.i => Symbolics.solve_for(J1_eqs[2], scsim.J1.sys.i)))

eqs15= substitute(eqs14, Dict(scsim.loop3.sys.iₘ => Symbolics.solve_for(loop3_eqs[1], scsim.loop3.sys.iₘ)))

eqs16 = substitute(eqs15, Dict(scsim.R1.sys.i => Symbolics.solve_for(R1_eqs[1], scsim.R1.sys.i)))

eqs17 = substitute(eqs16, Dict(dRθ => Symbolics.solve_for(obs[4], dRθ)))

eqs18 = substitute(eqs17, Dict(scsim.loop4.sys.iₘ => Symbolics.solve_for(loop4_eqs[1], scsim.loop4.sys.iₘ)))

eqs19 = Symbolics.simplify(eqs18)

sys = [eqs10, eqs19]

ddJ1 = get_variables(obs[1])[1]

ddC1 = get_variables(obs[2])[1]

dJ1 = states(model)[2]
dC1 = states(model)[4]

using HarmonicBalance

sys2 = substitute(sys, Dict(ddJ1 => d(scsim.J1.sys.θ,t,2)))

sys3 = substitute(sys2, Dict(ddC1 => d(scsim.C1.sys.θ,t,2)))

sys4 = substitute(sys3, Dict([dJ1 => d(scsim.J1.sys.θ,t), dC1 => d(scsim.C1.sys.θ,t)]))

diff_eq1 = DifferentialEquation(sys4, [scsim.J1.sys.θ, scsim.C1.sys.θ])

add_harmonic!(diff_eq1, scsim.J1.sys.θ, scsim.loop4.sys.ω)
add_harmonic!(diff_eq1, scsim.C1.sys.θ, scsim.loop4.sys.ω)

harmonic_eq = get_harmonic_equations(diff_eq1)

varied = scsim.loop4.sys.ω => 2*pi*100.0e+6*LinRange(0.0, 2.0, 200)

fixed = (scsim.loop4.sys.I => 3.0*I₀,
scsim.J1.sys.I0 => I₀,
scsim.J1.sys.R => R₀,
scsim.J1.sys.C => 0.01/βc,
scsim.J1.sys.L => 5.0/βL,
scsim.R1.sys.R => 50.0,
scsim.C1.sys.C => 2.0/βc,
scsim.L1.sys.L => 2.0/βL,
scsim.L2.sys.L => 100.0/βL,
scsim.M12.sys.L => 8.0/βL,
scsim.loop1.sys.Φₑ => 0.25*Φ₀,
scsim.loop2.sys.Φₑ => 0.25)

result = get_steady_states(harmonic_eq, varied, fixed)

plot(result, "sqrt(u1^2 + v1^2)")

sys4