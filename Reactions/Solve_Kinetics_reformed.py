from numpy import dot, log, zeros, array, exp
from Parameters import R


def solve_kinetics(Cl, Tl):

    Rf = zeros((2, 1))
    Rr = zeros((2, 1))

    kf1 = 3.951e18 * exp(-6864 / Tl)
    kf2 = 2.11e1

    C3 = 2.8898, - 3635.09, 0, 0
    C4 = 2.1211, - 8189.38, 0, - 0.007484
    C5 = 231.465, - 12092.10, - 36.7816, 0

    def f_K(C, T):
        C1, C2, C3, C4 = C
        return exp(C1 + C2 / T + C3 * log(T) + C4 * T)

    K3 = f_K(C3, Tl)
    K4 = f_K(C4, Tl)
    K5 = f_K(C5, Tl)

    K1 = K5/(K4*K3)
    K2 = K5/K4

    # actual concentration equilibrium
    Kee1 = (Cl[5] * Cl[6] / (Cl[0] * Cl[1] ** 2))
    Kee2 = (Cl[5] * Cl[7] / (Cl[0] * Cl[1] * Cl[2]))

    # Reverse reaction rate
    kr1 = kf1 / K1
    kr2 = kf2 / K2

    # forward and reverse reaction rates
    Rf1 = Cl[0] * (Cl[1] ** 2)  * kf1
    Rr1 = Cl[5] * Cl[6]         * kr1

    Rf2 = Cl[0] * Cl[1] * Cl[2] * kf2
    Rr2 = Cl[5] * Cl[7]         * kr2

    Rf[0, 0] = Rf1
    Rf[1, 0] = Rf2
    Rr[0, 0] = Rr1
    Rr[1, 0] = Rr2

    Ra = Rf.T - Rr.T  # overall reaction rate

    print(Ra)

    # actual eq constant minus numerical value at eq are equal
    Kt1 = Kee1 - K1
    Kt2 = Kee2 - K2

    # stoichiometric matrix
    ST = array([[-1, -2, 0, 0, 0, 1, 1, 0],
                [-1, -1, -1, 0, 0, 1, 0, 1]])

    Rgen = dot(ST.T, Ra.T).T[0]  # rate of generation vector

    # print(Rgen)

    return Rgen
