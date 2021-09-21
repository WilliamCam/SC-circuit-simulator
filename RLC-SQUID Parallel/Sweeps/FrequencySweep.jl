#Problem setup in order to vary drive frequency ωin

using UnPack

#Custom structs to keep number of mutable parameters low
struct df_sweep_consts
    Ib::Float64
    Iin::Float64
    Φe::Float64
end
mutable struct df_sweep_vars
    ω::Float64
    i::Vector{Float64}
end

struct df_sweep_params
    c::df_sweep_consts
    v::df_sweep_vars
end
"""
    RLC_para_fswp!(du,u,p,t)

Defines RLC-SQUID w/LPF circuit problem with variable drive frequency ωin
"""
function RLC_para_fswp!(du,u,p,t)
    @unpack c,v = p
    @unpack Ib, Iin, Φe = c
    @unpack ω, i = v

    i_mesh!(i, Ib, Iin, Φe, ω, t, u)




     du[1] = EoM_junction1(i, Gλ)
     du[2] = u[7]
     du[3] = u[8]
     du[4] = u[9]
     du[5] = EoM_junction5(i,G1)
     du[6] = u[10]

     du[7] = EoM_junction2(i, Cλ)
     du[8] = EoM_junction3(i, Cj1, I₀₁, Gj1, u)
     du[9] = EoM_junction4(i, Cj2, I₀₂, Gj2, u)
     du[10] = u[11]

     du[11] = EoM_junction6(i, C1, C2, G2, u)

end #end RLC_para_fswp
#initialising the junction current cache vector
junc_i = zeros(6)
u0 = zeros(11)

"""
    RLC_para_fswp_S21!(du,u,p,t)

Defines RLC-SQUID w/LPF circuit problem with variable drive frequency ωin.
solution vector u[] is extended to contain the circuit output current so that S21 parameters
can be determined
"""
function RLC_para_fswp_S21!(du,u,p,t)
    @unpack c,v = p
    @unpack Ib, Iin, Φe = c
    @unpack ω, i = v

    i_mesh!(i, Ib, Iin, Φe, ω, t, u) #determine currents flowing through each junction




     du[1] = EoM_junction1(i, Gλ) #dθ₁/dt
     du[2] = u[7] #dθ₂/dt
     du[3] = u[8] #dθ₃/dt
     du[4] = u[9] #dθ₄/dt
     du[5] = EoM_junction5(i,G1) #dθ₅/dt
     du[6] = u[10] #dθ₆/dt

     du[7] = EoM_junction2(i, Cλ) #d²θ₂/dt²
     du[8] = EoM_junction3(i, Cj1, I₀₁, Gj1, u) #d²θ₃/dt²
     du[9] = EoM_junction4(i, Cj2, I₀₂, Gj2, u) #d²θ₄/dt²
     du[10] = u[11] #d²θ₆/dt²

     du[11] = EoM_junction6(i, C1, C2, G2, u) #d³θ₆/dt³

end #end RLC_para_fswp_S21

u0_S21 = zeros(12)
