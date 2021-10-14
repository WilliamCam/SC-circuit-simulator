
using UnPack

#Custom structs to keep number of mutable parameters low
struct flux_sweep_consts
    Ib::Float64
    Iin::Float64
    ω::Float64
end
mutable struct flux_sweep_vars
    Φe::Float64
    i::Vector{Float64}
end

struct flux_sweep_params
    c::flux_sweep_consts
    v::flux_sweep_vars
end

"""
    RLC_para_fluxswp!(du,u,p,t)

Defines RLC-SQUID w/LPF circuit problem with variable flux bias Φₑ
"""
function RLC_para_v2_fluxswp!(du,u,p,t)
    @unpack c,v = p
    @unpack Ib, Iin, ω = c
    @unpack Φe, i = v

    i_mesh!(i, Ib, Iin, Φe, ω, t, u)




     du[1] = u[6]
     du[2] = u[7]
     du[3] = u[8]
     du[4] = u[9]
     du[5] = u[10]

     du[6] = d²θC₀(i,u,Ginp,C0)
     du[7] = d²θCλ(i,u,Cλ)
     du[8] = d²θ₃(i,u,Gj1,Cj1,I₀₁)
     du[9] = d²θ₄(i,u,Gj2,Cj2,I₀₂)
     du[10] = u[11]

     du[11] = d³θC₂(i,u,G1,G2,C1,C2)

end #end RLC_para_v2_fluxswp


#Custom structs to keep number of mutable parameters low
struct freq_sweep_consts
    Ib::Float64
    Iin::Float64
    Φe::Float64
end
mutable struct freq_sweep_vars
    ω::Float64
    i::Vector{Float64}
end

struct freq_sweep_params
    c::freq_sweep_consts
    v::freq_sweep_vars
end

function RLC_para_v2_fswp!(du,u,p,t)
    @unpack c,v = p
    @unpack Ib, Iin, Φe = c
    @unpack ω, i = v

    i_mesh!(i, Ib, Iin, Φe, ω, t, u)




     du[1] = u[6]
     du[2] = u[7]
     du[3] = u[8]
     du[4] = u[9]
     du[5] = u[10]

     du[6] = d²θC₀(i,u,Ginp,C0)
     du[7] = d²θCλ(i,u,Cλ)
     du[8] = d²θ₃(i,u,Gj1,Cj1,I₀₁)
     du[9] = d²θ₄(i,u,Gj2,Cj2,I₀₂)
     du[10] = u[11]

     du[11] = d³θC₂(i,u,G1,G2,C1,C2)

end #end RLC_para_v2_fswp

function RLC_para_v2_fluxswp_TEST!(du,u,p,t)
    @unpack c,v = p
    @unpack Ib, Iin, ω = c
    @unpack Φe, i = v

    i_mesh!(i, Ib, Iin, Φe, ω, t, u)

     du[1] = u[6]
     du[2] = u[7]
     du[3] = u[8]
     du[4] = u[9]
     du[5] = u[10]

     du[6] = d²θC₀(i,u,Ginp,C0)
     du[7] = d²θCλ(i,u,Cλ)
     du[8] = d²θ₃(i,u,Gj1,Cj1,I₀₁)
     du[9] = d²θ₄(i,u,Gj2,Cj2,I₀₂)
     du[10] = u[11]

     du[11] = d³θC₂(i,u,G1,G2,C1,C2)
     u[12] = i[1]
     u[13] = i[2]

end #end RLC_para_v2_fswp
