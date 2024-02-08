import numpy as np


def solve_driving_force(x, y, Tl, a0, alpha, KH, H_CO2_mix, P):

    # IDAES Parameters for Psat H2O
    f_Psat_H2O = lambda T: np.exp(72.55 + -7206.70 / T + -7.1385 * np.log(T) + 4.05e-6 * T ** 2)

    # From Gabrielsen: A Model for Estimating CO2 Solubility in Aqueous Alkanolamines
    K_CO2 = np.exp(30.96 + -10584 / Tl + -7.187 * a0 * alpha)

    Psat_H2O = f_Psat_H2O(Tl)

    Pv_CO2 = y[0] * P
    # From Gabrielsen: A Model for Estimating CO2 Solubility in Aqueous Alkanolamines
    Pl_CO2 = K_CO2 * x[0] * a0 * alpha / (a0 * (1 - 2 * alpha)) ** 2

    Pv_H2O = y[1] * P*.9
    Pl_H2O = x[2]*Psat_H2O

    DF_CO2 = (Pv_CO2 - Pl_CO2)*KH
    DF_H2O = (Pv_H2O - Pl_H2O)

    return DF_CO2, DF_H2O, [DF_CO2, Pv_CO2, Pl_CO2, H_CO2_mix, DF_H2O, Pv_H2O, Pl_H2O, Psat_H2O]
