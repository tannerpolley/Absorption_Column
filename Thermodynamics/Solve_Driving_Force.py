import numpy as np
from Parameters import f_interp


def solve_driving_force(x, y, Tl, a0, alpha, KH, H_CO2_mix, P, zi):

    # Cli, Cvi, K, Z, Φ = solve_vle(Cl, Cv, Tl, Tv, df_param)

    # IDAES Parameters for Psat H2O
    f_Psat_H2O = lambda T: np.exp(72.55 + -7206.70 / T + -7.1385 * np.log(T) + 4.05e-6 * T ** 2)
    K_CO2 = np.exp(30.96 + -10584 / Tl + -7.187 * a0 * alpha)

    Psat_H2O = f_Psat_H2O(Tl)
    DF_H2O_IDAES = f_interp('DF_H2O', zi)
    Pv_H2O_IDAES = f_interp('Pv_H2O', zi)
    Pl_H2O_IDAES = f_interp('Pl_H2O', zi)
    Psat_H2O_IDAES = f_interp('Psat_H2O', zi)
    x_H2O_true = f_interp('x_H2O_true', zi)

    Pv_CO2 = y[0] * P

    Pl_CO2 = K_CO2 * x[0] * a0 * alpha / (a0 * (1 - 2 * alpha)) ** 2
    # Pl_CO2 = H_CO2_mix * Cl_true[0]

    Pv_H2O = y[1] * P
    # Pv_H2O = Pv_H2O_IDAES
    Pl_H2O = x[2]*Psat_H2O

    DF_CO2 = (Pv_CO2 - Pl_CO2)*KH
    DF_H2O = (Pv_H2O - Pl_H2O)

    # if zi < 1:
    #     DF_H2O = DF_H2O/2

    return DF_CO2, DF_H2O, [DF_CO2, Pv_CO2, Pl_CO2, H_CO2_mix, DF_H2O, Pv_H2O, Pl_H2O, Psat_H2O, DF_H2O_IDAES, Pv_H2O_IDAES, Pl_H2O_IDAES, Psat_H2O_IDAES]
