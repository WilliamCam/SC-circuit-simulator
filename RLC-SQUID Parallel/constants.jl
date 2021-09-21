#Constants

#some physical constants
const phi0 = 2.067833848e-15
const pi = 3.14159265359

#SQUID parameters

const I₀₁ = 10.0e-6 #junction critical current
const I₀₂ = 10.0e-6

const Ls = 100.0e-12 #SQUID loop inductance
const Gj1 = 0.1 #junction admittance
const Gj2 = 0.1

const βc₁ = 0.1 #Stewart-McCumber parameter (defines junction capacitance)
const βc₂ = 0.1

const Cj1 = βc₁*phi0/(2*pi)*Gj1^2/(I₀₁) #junction capacitance
const Cj2 = βc₂*phi0/(2*pi)*Gj2^2/(I₀₂)

const ωc = 2*pi*I₀₁/(phi0*Gj1) #characteristic Josephsons frequency

#Input coil parameters

const Linp = 500.0e-9 #input inductance
const Minp = phi0/0.525e-6 #input coil coupling

#resonator parameters, start with small Q to reduce simulation time.
const ωλ = 200.0e+6*2*pi #oscillation frequency
const Qλ = 1e6 #Quality factor
const Gλ = 1.0/2.9 #motional admittance
const Rλ = 1/Gλ #motional resistance
const Cλ = Gλ/(Qλ*ωλ)
const Lλ = Qλ/(Gλ*ωλ)

#low pass filter parameters
const G1 = 1.0e-3
const G2 = 1.0e-4
const C1 = 0.32e-12
const C2 = 0.032e-12;

#oscillation parameters
const tau = Qλ/ωλ
