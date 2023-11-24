import numpy as np
from Parameters import vapor_species


def liquid_diffusivity(Tl, C_MEA, df_param):

    C_MEA = C_MEA * 1e-3

    # Get Diffusivity of Liquid
    a, b, c, d, e = df_param['Dl_CO2'].dropna().to_numpy()

    Dl_CO2 = (a + b * C_MEA + c * C_MEA**2) * np.exp((d + (e * C_MEA))/Tl)

    return Dl_CO2


def vapor_diffusivity(Tv, y, P,  df_param):


    coefficients = df_param['Vapor Diffusivity'].dropna().to_numpy()

    Dv = (coefficients * Tv ** 1.75) / P

    Dv_T = 0
    for i in range(len(vapor_species)):
        Dv_T += y[i] * Dv[i]

    Dv_CO2, Dv_H2O, Dv_N2, Dv_O2 = Dv

    return Dv_CO2, Dv_H2O, Dv_N2, Dv_O2, Dv_T
