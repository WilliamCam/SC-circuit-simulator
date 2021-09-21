#Problem setup in order to vary flux bias Φₑ

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
function RLC_para_fluxswp!(du,u,p,t)
    @unpack c,v = p
    @unpack Ib, Iin, ω = c
    @unpack Φe, i = v

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
