#this file contains the governing circuit equations of motion for the phase difference θ and currnet i of each component
# for the RLC-SQUID Parallel circuit configuration with low pass filter v2
"""
    i_mesh!(i, Ib, Iin, Φe, ωin, t, u)

computes current at time t flowing through each junction.
Result is stored in cache vector i.
"""
function i_mesh!(i, Ib, Iin, Φe, ωin, t, u)
    i[1] = (2*G1*G2*Gλ*Lλ*pi*(Iin*sin(ωin*t)*Linp*Ls - 2*Iin*sin(ωin*t)*Minp^2 + 2*Minp*Φe) +
         phi0*(-(G1*G2*Gλ*(Ls*(Linp + Lλ) - 2*Minp^2)*u[1]) - 2*G1*G2*Gλ*Lλ*Minp*u[3] +
            2*G1*G2*Gλ*Lλ*Minp*u[5] + G1*G2*Gλ*Linp*Ls*u[2] - 2*G1*G2*Gλ*Minp^2*u[2] +
            4*C2*G1*Gλ*Lλ*Minp*u[10] + 2*C1*G2*Gλ*Lλ*Minp*u[10] +
            Cλ*G1*G2*Linp*Ls*u[7] - 2*Cλ*G1*G2*Minp^2*u[7] +
            2*C1*C2*Gλ*Lλ*Minp*u[11]))/(2*G1*G2*Gλ*Lλ*(Linp*Ls - 2*Minp^2)*pi)

    i[2] = (phi0*(Gλ*u[1] - Gλ*u[2] - Cλ*u[7]))/(2*Gλ*Lλ*pi)

    i[3] =   (2*G1*G2*Linp*pi*Φe - G1*G2*Minp*phi0*u[1] +
         Linp*phi0*(-(G1*G2*u[3]) + G1*G2*u[5] + (2*C2*G1 + C1*G2)*u[10] +
         C1*C2*u[11]))/(G1*G2*(Linp*Ls - 2*Minp^2)*pi)

    i[4] = (phi0*(-(G1*G2*u[4]) + G1*G2*u[5] + (2*C2*G1 + C1*G2)*u[10] + C1*C2*u[11]))/(G1*G2*Ls*pi)

    i[5] =  (G1*G2*Ls*pi*(Ib*Linp*Ls - 2*Ib*Minp^2 - 2*Linp*Φe) + G1*G2*Ls*Minp*phi0*u[1] +
        G1*G2*Linp*Ls*phi0*u[3] + G1*G2*(Linp*Ls - 2*Minp^2)*phi0*u[4] -
        2*(Linp*Ls - Minp^2)*phi0*(G1*G2*u[5] + (2*C2*G1 + C1*G2)*u[10] +
        C1*C2*u[11]))/(G1*G2*Ls*(Linp*Ls - 2*Minp^2)*pi)

end #end i_mesh!

#equations of motion for phase across each junction

function d²θC₀(i,u,G,C)
    2*pi/(phi0*C)*(i[1]-phi0*G/(2*pi)*u[6])
end #end d²θC₀

function d²θCλ(i,u,C)
    i[2]*2*pi/(phi0*C)
end #end d²θCλ

function d²θ₃(i,u,G,C,I0)
    2*pi/(phi0*C)*(i[3] - I0*sin(u[3]) - G*phi0*u[8]/(2*pi))
end #end d²θ₃

function d²θ₄(i,u,G,C,I0)
    2*pi/(phi0*C)*(i[4] - I0*sin(u[4]) - G*phi0*u[9]/(2*pi))
end #end d²θ₃

function d³θC₂(i,u,G1,G2,C1,C2)
    G2/(C1*C2)*(2*pi*i[5]/phi0-C1*u[11]-C2*u[11])
end #end d³θC₂
