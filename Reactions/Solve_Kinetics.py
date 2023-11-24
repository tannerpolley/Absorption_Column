from numpy import log, exp
import numpy as np


def solve_kinetics(Cl, Tl):

    kf1 = 3.1732e9 * exp(-4936.6 / Tl)
    kf2 = 1.51e5 * exp(-3092.9 / Tl)


    # C1 = 233.4, -3410, -36.8, 0
    # C2 = 176.72, -2909, -28.46, 0
    C1 = 164.039636, -707.0056712, -26.40136817, 0
    C2 = 366.061867998774, -13326.25411, -55.68643292, 0

    def f_logK(C, T):
        C1, C2, C3, C4 = C
        return C1 + C2 / T + C3 * log(T) + C4 * T + log(1e-3)

    K1 = exp(f_logK(C1, Tl))
    K2 = exp(f_logK(C2, Tl))

    # Reverse rate constant
    kr1 = kf1 / K1
    kr2 = kf2 / K2

    # forward and reverse reaction rates

    Rf1 = Cl[0] * (Cl[1] ** 2) * kf1
    Rr1 = Cl[5] * Cl[6] * kr1

    Rf2 = Cl[0] * Cl[1] * Cl[2] * kf2
    Rr2 = Cl[5] * Cl[7] * kr2

    R1 = Rf1 - Rr1
    R2 = Rf2 - Rr2

    ST_1 = [-1, -2, 0, 0, 0, 1, 1, 0]  # 1
    ST_2 = [-1, -1, -1, 0, 0, 1, 0, 1]  # 2

    R_gen = np.zeros(len(Cl))

    for i in range(len(Cl)):
        R_gen[i] = (ST_1[i]*R1 + ST_2[i]*R2)

    return R_gen
