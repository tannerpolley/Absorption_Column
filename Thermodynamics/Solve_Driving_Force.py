import numpy as np
from Parameters import f_DF_H2O


def solve_driving_force(Cl_true, x_true, y, Tl, Tv, KH, H_CO2_mix, P, zi):

    # Cli, Cvi, K, Z, Φ = solve_vle(Cl, Cv, Tl, Tv, df_param)

    f_Psat_H2O = lambda T: np.exp(73.649 + -7258.2 / T + -7.3037 * np.log(T) + 4.1653e-6 * T ** 2)

    Psat_H2O = f_Psat_H2O(Tl)

    Pv_CO2 = y[0] * P
    Pl_CO2 = H_CO2_mix * Cl_true[0]

    Pv_H2O = y[1] * P
    Pl_H2O = x_true[2]*Psat_H2O

    DF_CO2 = (Pv_CO2 - Pl_CO2)*KH
    DF_H2O = (Pv_H2O - Pl_H2O)

    DF_H2O_f = f_DF_H2O(zi)

    DF_CO2 = DF_CO2

    if zi == 0:
        DF_H2O = DF_H2O_f
    else:
        DF_H2O = DF_H2O

    return DF_CO2, DF_H2O_f, [DF_CO2, DF_H2O, Pv_CO2, Pl_CO2, Pv_H2O, Pl_H2O, H_CO2_mix, Psat_H2O]
