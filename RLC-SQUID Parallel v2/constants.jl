#This file contains all constants used in the RLC-SQUID Parallel v2 simulation

#some physical constants
const phi0 = 2.067833848e-15
const pi = 3.14159265359

#SQUID parameters

const I₀₁ = 10.0e-6 #junction critical current
const I₀₂ = 10.0e-6

const Ls = 100.0e-12 #SQUID loop inductance
const Gj1 = 1/10 #junction admittance
const Gj2 = 1/10

const βc₁ = 0.1 #Stewart-McCumber parameter (defines junction capacitance)
const βc₂ = 0.1

const Cj1 = βc₁*phi0/(2*pi)*Gj1^2/(I₀₁) #junction capacitance
const Cj2 = βc₂*phi0/(2*pi)*Gj2^2/(I₀₂)

const ωc = 2*pi*I₀₁/(phi0*Gj1) #characteristic Josephsons frequency

#Input coil parameters
const Ginp = 1/500
const Linp = 400.0e-9 #input inductance
const Minp = phi0/0.525e-6 #input coil coupling
const βL = 2*I₀₁*Ls/phi0

#resonator parameters, start with small Q to reduce simulation time.
const ωλ = 100.0e+6*2*pi
const Qλ = 1e3 #Quality factor
const Gλ = 1.0/2.9 #motional admittance
const Rλ = 1/Gλ #motional resistance
const Cλ = Gλ/(Qλ*ωλ) #motional capacitacne (determined by ω and Q)
const Lλ = Qλ/(Gλ*ωλ) #motional inductance (determined by ω and Q)
const C0 = 1.0e-12#shunt capacitance of electrodes

#low pass filter parameters ### Large G impacts performance
const G1 = 1.0e-2
const G2 = 1.0e-3
const C1 = 0.32e-11
const C2 = 0.032e-11;

#oscillation parameters
const τ = 2*Qλ/ωλ
