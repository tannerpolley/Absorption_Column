import numpy as np
exp = np.exp
array = np.array
sum = np.sum
π = np.pi
log = np.log
import numdifftools as nd

k = 1.3806e-23 # J/K
N_A = 6.0221e23 # 1/mol
R = 8.314 # J/mol-K
η = .402 # Packing fraction
v = 6.3e-5 # m^3/mol

a_ni = np.array([[0.910563,-0.308402,-0.090615],
                 [0.636128,0.186053,0.452784],
                 [2.686135,-2.503005,0.596270],
                 [-26.547362,21.419794,-1.724183],
                 [97.759209,-65.255885,-4.130211],
                 [-159.591541,83.318680,13.776632],
                 [91.297774,-33.746923,-8.672847]])
a_ni = a_ni.T

b_ni = np.array([[0.724095,-0.575550,0.097688],
                 [2.238279,0.699510,-0.255757],
                 [-4.002585,3.892567,-9.155856],
                 [-21.003577,-17.215472,20.642076],
                 [26.855641,192.672264,-38.804430],
                 [206.551338,-161.826462,93.626774],
                 [-355.602356,-165.207693,-29.666906]])
b_ni = b_ni.T

T = 315 # Temperature (K)

# x = array([.05, .3, .65]) # Mole fraction
# m = array([2.0729, 3.0353, 1.9599]) # Number of segments
# σ = array([2.7852, 3.0435, 2.362]) # Temperature-Independent segment diameter σ_i (Aᵒ)
# ϵ_k = array([169.21, 277.174, 279.42]) # Depth of pair potential / Boltzmann constant (K)
# k_ij = array([[0.00E+00,.16,.065],
#               [.16,0.00E+00,-.18],
#               [.065,-.18,0.00E+00]])
# κ_AB = array([0, .037470, .2039])
# ϵ_AB_k = array([0, 2586.3, 2059.28])


T = 233.15 # Temperature (K)
x = array([.1, .3, .6]) # Mole fraction
m = array([1, 1.6069, 2.0020]) # Number of segments
σ = array([3.7039, 3.5206, 3.6184]) # Temperature-Independent segment diameter σ_i (Aᵒ)
ϵ_k = array([150.03, 191.42, 208.11]) # Depth of pair potential / Boltzmann constant (K)
k_ij = array([[0.00E+00,3.00E-04,1.15E-02],
              [3.00E-04,0.00E+00,5.10E-03],
              [1.15E-02,5.10E-03,0.00E+00]])

κ_AB = array([0, 0, 0])
ϵ_AB_k = array([0, 0, 0])

def f_a_res(η):

    d = σ*(1-.12*exp(-3*ϵ_k/T))

    ρ = 6/π*η*(sum([x[i]*m[i]*d[i]**3 for i in range(3)]))**(-1)

    ξ = array([])
    for n in range(4):
        Σ = 0
        for i in range(3):
            Σ += x[i]*m[i]*d[i]**n
        ξ = np.append(ξ, π/6*ρ*Σ)

    g_hs_ij = np.zeros((3, 3))
    for i in range(3):
        for j in range(3):
            g_hs_ij[i][j] = 1/(1 - ξ[3]) + (d[i]*d[j]/(d[i]+d[j]))*3*ξ[2]/(1-ξ[3])**2 + (d[i]*d[j]/(d[i]+d[j]))**2 * 2*ξ[2]**2/(1-ξ[3])**3

    ã_hs = 1/ξ[0]*(3*ξ[1]*ξ[2]/(1-ξ[3]) + ξ[2]**3/(ξ[3]*(1-ξ[3])**2) + (ξ[2]**3/ξ[3]**2 - ξ[0])*log(1-ξ[3]))

    m̄ = sum(x*m)

    ã_hc = m̄*ã_hs - sum([x[i]*(m[i] - 1)*log(g_hs_ij[i][i]) for i in range(3)])

    a = a_ni[0]+ (m̄-1)/m̄*a_ni[1] + (m̄-1)/m̄*(m̄-2)/m̄*a_ni[2]
    b = b_ni[0]+ (m̄-1)/m̄*b_ni[1] + (m̄-1)/m̄*(m̄-2)/m̄*b_ni[2]

    I1 = sum([a[i]*η**i for i in range(7)])
    I2 = sum([b[i]*η**i for i in range(7)])

    σ_ij = np.zeros((3, 3))
    for i in range(3):
        for j in range(3):
            σ_ij[i][j] = 1 / 2 * (σ[i] + σ[j])

    ϵ_ij = np.zeros((3, 3))
    for i in range(3):
        for j in range(3):
            ϵ_ij[i][j] = (ϵ_k[i] * ϵ_k[j]) ** (1 / 2) * (1 - k_ij[i][j])

    Σ_i = 0
    for i in range(3):
        Σ_j = 0
        for j in range(3):
            Σ_j += x[i]*x[j]*m[i]*m[j]*(ϵ_ij[i][j]/T)*σ_ij[i][j]**3
        Σ_i += Σ_j
    Σ_1 = Σ_i

    Σ_i = 0
    for i in range(3):
        Σ_j = 0
        for j in range(3):
            Σ_j += x[i]*x[j]*m[i]*m[j]*(ϵ_ij[i][j]/T)**2*σ_ij[i][j]**3
        Σ_i += Σ_j
    Σ_2 = Σ_i

    C1 = (1 + m̄*(8*η - 2*η**2)/(1-η)**4 + (1-m̄)*(20*η - 27*η**2 + 12*η**3 - 2*η**4)/((1-η)*(2-η))**2)**-1

    ã_disp = -2*π*ρ*I1*Σ_1 - π*ρ*m̄*C1*I2*Σ_2

    κ_AB_ij = np.zeros((3, 3))
    for i in range(3):
        for j in range(3):
            κ_AB_ij[i][j] = (κ_AB[i] * κ_AB[j]) ** (1 / 2) * ((σ[i] * σ[j]) / (1 / 2 * (σ[i] * σ[j]))) ** 3

    ϵ_AB_ij = np.zeros((3, 3))
    for i in range(3):
        for j in range(3):
            ϵ_AB_ij[i][j] = (ϵ_AB_k[1] + ϵ_AB_k[2]) / 2

    d_ij = np.zeros((3, 3))
    for i in range(3):
        for j in range(3):
            d_ij[i][j] = 1/2*(d[i] + d[j])

    Δ_AB_ij = np.zeros((3, 3))
    for i in range(3):
        for j in range(3):
            Δ_AB_ij[i][j] = d_ij[i][j]**3*g_hs_ij[i][j]*κ_AB_ij[i][j]*(exp(ϵ_AB_ij[i][j]/T) - 1)


    def XA_find(XA_guess, n, Δ_AB_ij, ρ, x):
        m = int(XA_guess.shape[1] / n)
        X_Ai = np.zeros_like(XA_guess)
        AB_matrix = np.asarray([[0., 1.],
                                [1., 0.]])
        Σ_2 = np.zeros((n,), dtype='float_')
        XA = np.zeros_like(XA_guess)

        for i in range(n):
            Σ_2 = 0 * Σ_2
            for j in range(n):
                Σ_2 += ρ * x[j] * (XA_guess[j, :] @ (Δ_AB_ij[i][j] * AB_matrix))
            XA[i, :] = 1 / (1 + Σ_2)

        return XA


    ncA = 2
    a_sites = 2
    iA = [1, 2]
    XA = np.zeros((ncA, a_sites), dtype='float_')

    ctr = 0
    dif = 1000.
    XA_old = np.copy(XA)
    while (ctr < 500) and (dif > 1e-9):
        ctr += 1
        XA = XA_find(XA, ncA, Δ_AB_ij, ρ, x[iA])
        dif = np.sum(abs(XA - XA_old))
        XA_old[:] = XA
    XA = XA.flatten()

    ã_assoc = 0
    for i in range(3):
        Σ_1 = 0
        for j in range(len(XA)):
            Σ_1 += log(XA[j] - 1 / 2 * XA[j] + 1 / 2)

        ã_assoc += x[i] * Σ_1

    ã_res = ã_hc + ã_disp + ã_assoc
    # print(ã_hc, ã_disp)
    return ã_res


da_dn = nd.Derivative(f_a_res)

# print(da_dn(.402))
Z = 1 + .402*da_dn(.402)
print(Z)