#Equations of motion for RLC-SQUID coupled open loop circuit with 2nd order low pass filter

"""
    i_mesh!(i, Ib, Iin, Φe, ωin, t, u)

computes current at time t flowing through each junction.
Result is stored in cache vector i.
"""
function i_mesh!(i, Ib, Iin, Φe, ωin, t, u)
     i[1] = (2.0*Iin*sin(ωin*t)*Linp*Ls*pi - 4.0*Iin*sin(ωin*t)*Minp^2.0*pi + 4.0*Minp*pi*Φe - Ls*phi0*u[1] - Ls*phi0*u[2] - 2.0*Minp*phi0*u[3] +
        2.0*Minp*phi0*u[5] + 2.0*Minp*phi0*(C2*Rλ*u[10]+u[6]))/(2.0*(Linp*Ls + Ls*Lλ - 2.0*Minp^2)*pi)

    i[2] = i[1]

    i[3] = (-2.0*Iin*sin(ωin*t)*Lλ*Minp*pi + 2.0*Linp*pi*Φe + 2.0*Lλ*pi*Φe - Minp*phi0*u[1] - Minp*phi0*u[2] - Linp*phi0*u[3] -
        Lλ*phi0*u[3] + Linp*phi0*u[5] + Lλ*phi0*u[5] + Linp*phi0*(C2*Rλ*u[10]+u[6]) + Lλ*phi0*(C2*Rλ*u[10]+u[6]))/((Linp*Ls + Ls*Lλ - 2.0*Minp^2.0)*pi)

    i[4] = -1.0*((phi0*(u[4] - u[5] - (C2*Rλ*u[10]+u[6])))/(Ls*pi))

    i[5] = (1.0/(Ls*(Linp*Ls + Ls*Lλ - 2.0*Minp^2.0)*pi))*(Ib*Linp*Ls^2.0*pi + Ib*Ls^2.0*Lλ*pi + 2.0*Iin*sin(ωin*t)*Ls*Lλ*Minp*pi
        - 2.0*Ib*Ls*Minp^2.0*pi - 2.0*Linp*Ls*pi*Φe - 2.0*Ls*Lλ*pi*Φe + Ls*Minp*phi0*u[1] + Ls*Minp*phi0*u[2]
        + Linp*Ls*phi0*u[3] + Ls*Lλ*phi0*u[3] + Linp*Ls*phi0*u[4] + Ls*Lλ*phi0*u[4] - 2.0*Minp^2.0*phi0*u[4]
        - 2.0*Linp*Ls*phi0*u[5] - 2.0*Ls*Lλ*phi0*u[5] + 2.0*Minp^2*phi0*u[5] - 2.0*Linp*Ls*phi0*(C2*Rλ*u[10]+u[6]) - 2.0*Ls*Lλ*phi0*(C2*Rλ*u[10]+u[6])
        + 2.0*Minp^2*phi0*(C2*Rλ*u[10]+u[6]))

    i[6] = i[5]
end


"""
    EoM_junction1(i,G)

Equation of motion for of junction 1
"""
function EoM_junction1(i, G)
    2*pi*i[1]/(phi0*G)
end
"""
    EoM_junction1(i,G)

Equation of motion for of junction 2
"""
function EoM_junction2(i,C)
     2*pi/(phi0*C)*i[2]
end

"""
    EoM_junction1(i,G)

Equation of motion for of junction 3
"""
function EoM_junction3(i,C,I0,G,u)
   2*pi/(phi0*C)*(i[3] - I0*sin(u[3]) - G*phi0*u[8]/(2*pi))
end
"""
    EoM_junction1(i,G)

Equation of motion for of junction 4
"""
function EoM_junction4(i,C,I0,G,u)
   2*pi/(phi0*C)*(i[4] - I0*sin(u[4]) - G*phi0*u[9]/(2*pi))
end
"""
    EoM_junction1(i,G)

Equation of motion for of junction 5
"""
function EoM_junction5(i, G)
    2*pi*i[5]/(phi0*G)
end
"""
    EoM_junction1(i,G)

Equation of motion for of junction 6
"""
function EoM_junction6(i,C1,C2,G,u)
    (2.0*pi*G)/(phi0*C1*C2)*(i[6]-phi0*C1/(2.0*pi)*u[11]-phi0*C2/(2.0*pi)*u[11])
end
