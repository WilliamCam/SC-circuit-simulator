using DifferentialEquations, Plots, Statistics
include("C://Users//21958742//GitHub//SC-circuit-simulator//program//SCsim.jl")
scsim.open_file("C://Users//21958742//GitHub//SC-circuit-simulator//examples//LC-oscillator.jld2")

model, u0, old_model = scsim.build_circuit()

using HarmonicBalance, ModelingToolkit
eqs = equations(model)
sts = states(model)

#eqs = substitute(eqs, Dict(ps))

eqs1 = [scsim.D2(scsim.C1.sys.θ) ~ (6.283185307179586*(1 + scsim.C1.sys.k * (scsim.C1.sys.i^2))*scsim.C1.sys.i) / (2.067e-15*scsim.C1.sys.C)]
#eqs1 = substitute(eqs1, Dict(scsim.C1.sys.i => scsim.Rin.sys.i - scsim.loop2.sys.iₘ))
append!(eqs1,eqs[3:end])
eqs1 = substitute(eqs1, Dict(sts[1]=>scsim.D(scsim.C1.sys.θ)))

diff_eq = DifferentialEquation(eqs1[1:3], [scsim.C1.sys.θ, scsim.V1.sys.θ, scsim.C1.sys.i])

add_harmonic!(diff_eq, scsim.C1.sys.θ, scsim.V1.sys.ω)
add_harmonic!(diff_eq, scsim.V1.sys.θ, scsim.V1.sys.ω)
add_harmonic!(diff_eq, scsim.C1.sys.i, scsim.V1.sys.ω)
add_harmonic!(diff_eq, scsim.loop2.sys.iₘ, scsim.V1.sys.ω)

harmonic_eq = get_harmonic_equations(diff_eq)

varied = scsim.V1.sys.ω => LinRange(0.9e+3, 1.1e+3, 100)
fixed = (
    scsim.C1.sys.C => 1.0e-3,
    scsim.L1.sys.L => 1.0e-3,
    scsim.V1.sys.V => 1.0,
    scsim.L1.sys.k => 1.0,
    scsim.C1.sys.k => 1.0,
    scsim.loop2.sys.Φₑ=>1.0,
    scsim.loop2.sys.iₘ=>1.0,

)

solutions = get_steady_states(harmonic_eq, varied, fixed)