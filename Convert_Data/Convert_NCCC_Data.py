import pandas as pd
import numpy as np
from Parameters import R, Tv_in, P, A, rho_l, Tl_0
from Properties.Density import liquid_density


def convert_NCCC_data(X, df_param):

    m_T_l, m_T_v, alpha, w_MEA, y_CO2 = X

    species = ['CO2', 'MEA', 'H2O', 'N2', 'O2', 'T']
    columns_l = ['w', 'm', 'MW', 'n', 'u', 'V', 'C', 'x']
    columns_v = ['w', 'm', 'MW', 'n', 'u', 'V', 'C', 'y']

    # Molecular Weights
    MW_CO2 = 44.0095 / 1000  # kg/mol
    MW_MEA = 61.0831 / 1000  # kg/mol
    MW_H2O = 18.0153 / 1000  # kg/mol
    MW_N2 = 28.0135 / 1000  # kg/mol
    MW_O2 = 31.9988 / 1000  # kg/mol

    MWs_2 = np.array([MW_CO2, MW_MEA, MW_H2O, MW_N2, MW_O2, 0])

    alpha_O2_N2 = 0.08485753604
    alpha_H2O_N2 = 0.1010145181

    # Liquid Calculations

    rho_l = liquid_density(Tl_0, w_MEA, alpha, df_param)

    # Find Liquid Mass Flow Rates
    m_MEA_l = w_MEA * m_T_l  # kg/s
    m_CO2_l = m_MEA_l * alpha / MW_MEA * MW_CO2  # kg/s
    m_H2O_l = m_T_l - m_MEA_l - m_CO2_l  # kg/s
    m_N2_l = 0
    m_O2_l = 0

    # Find Liquid Molar Flow Rates
    n_CO2_l = m_CO2_l / MW_CO2  # mole/s
    n_MEA_l = m_MEA_l / MW_MEA  # mole/s
    n_H2O_l = m_H2O_l / MW_H2O  # mole/s
    n_N2_l = m_N2_l / MW_N2  # mole/s
    n_O2_l = m_O2_l / MW_N2  # mole/s
    n_T_l = n_MEA_l + n_CO2_l + n_H2O_l  # mole/s

    # Find Liquid Mole Fractions
    x_CO2 = n_CO2_l / n_T_l
    x_MEA = n_MEA_l / n_T_l
    x_H2O = n_H2O_l / n_T_l
    x_N2 = n_N2_l / n_T_l
    x_O2 = n_O2_l / n_T_l
    x_T = x_CO2 + x_MEA + x_H2O + x_N2 + x_O2


    # Find Liquid Weight Fractions
    w_CO2_l = m_CO2_l / m_T_l
    w_MEA_l = m_MEA_l / m_T_l
    w_H2O_l = m_H2O_l / m_T_l
    w_N2_l = m_N2_l / m_T_l
    w_O2_l = m_O2_l / m_T_l
    w_T_l = w_CO2_l + w_MEA_l + w_H2O_l + w_N2_l + w_O2_l

    # Find liquid Volumetric Flow Rate
    x_ls = [x_CO2, x_MEA, x_H2O]
    MW_l = 0
    for i in range(len(x_ls)):
        MW_l += x_ls[i] * MWs_2[:-1][i]  # kg/mol

    V_l = MW_l/rho_l

    # Find Liquid Velocity
    u_l = V_l / A * n_T_l

    # Find Liquid Concentrations
    C_CO2_l = x_CO2 / V_l
    C_MEA_l = x_MEA / V_l
    C_H2O_l = x_H2O / V_l
    C_N2_l = x_N2 / V_l
    C_O2_l = x_O2 / V_l
    C_T_l = C_CO2_l + C_MEA_l + C_H2O_l + C_N2_l + C_O2_l

    # Make arrays for DataFrame
    C_l = np.array([C_CO2_l, C_MEA_l, C_H2O_l, C_N2_l, C_O2_l, C_T_l])
    w_l = np.array([w_CO2_l, w_MEA_l, w_H2O_l, w_N2_l, w_O2_l, w_T_l])
    m_l = np.array([m_CO2_l, m_MEA_l, m_H2O_l, m_N2_l, m_O2_l, m_T_l])
    n_l = np.array([n_CO2_l, n_MEA_l, n_H2O_l, n_N2_l, n_O2_l, n_T_l])
    u_l = np.array([u_l, u_l, u_l, u_l, u_l, u_l])
    V_l = np.array([V_l, V_l, V_l, V_l, V_l, V_l])
    x = np.array([x_CO2, x_MEA, x_H2O, x_N2, x_O2, x_T])

    liquid_data = [w_l, m_l, MWs_2, n_l, u_l, V_l, C_l, x]

    liquid_dict = {}
    for i, c in enumerate(columns_l):
        liquid_dict[c] = liquid_data[i]

    liquid = pd.DataFrame(liquid_dict, index=species)
    liquid.style.set_caption("liquid")

    # Vapor Calculations

    # Find Vapor Mole Fractions
    y_N2 = (1 - y_CO2) / (1 + alpha_H2O_N2 + alpha_O2_N2)
    y_O2 = y_N2 * alpha_O2_N2
    y_H2O = y_N2 * alpha_H2O_N2
    y_MEA = 0
    y_T = y_CO2 + y_MEA + y_H2O + y_N2 + y_O2
    sigma = y_N2 * MW_N2 + y_O2 * MW_O2 + y_CO2 * MW_CO2 + y_H2O * MW_H2O

    # Find Vapor Mass Flow Rates

    w_CO2_v = y_CO2 * MW_CO2 / sigma
    w_MEA_v = y_MEA * MW_MEA / sigma
    w_H2O_v = y_H2O * MW_H2O / sigma
    w_N2_v = y_N2 * MW_N2 / sigma
    w_O2_v = y_O2 * MW_O2 / sigma

    m_CO2_v = w_CO2_v * m_T_v
    m_MEA_v = w_MEA_v * m_T_v
    m_H2O_v = w_H2O_v * m_T_v
    m_N2_v = w_N2_v * m_T_v
    m_O2_v = w_O2_v * m_T_v



    # Find Vapor Molar Flow Rates
    n_CO2_v = m_CO2_v / MW_CO2  # mole/s
    n_MEA_v = m_MEA_v / MW_MEA  # mole/s
    n_H2O_v = m_H2O_v / MW_H2O  # mole/s
    n_N2_v = m_N2_v / MW_N2  # mole/s
    n_O2_v = m_O2_v / MW_O2  # mole/s
    n_T_v = n_CO2_v + n_MEA_v + n_H2O_v + n_N2_v + n_O2_v  # mole/s

    # Find Liquid Weight Fractions
    w_CO2_v = m_CO2_v / m_T_v
    w_MEA_v = m_MEA_v / m_T_v
    w_H2O_v = m_H2O_v / m_T_v
    w_N2_v = m_N2_v / m_T_v
    w_O2_v = m_O2_v / m_T_v
    w_T_v = w_CO2_v + w_MEA_v + w_H2O_v + w_N2_v + w_O2_v

    # Find Vapor Volumetric Flow Rate
    V_v = n_T_v * R * Tv_in / P

    # Find Vapor Velocity
    u_v = V_v / A  # m/s

    # Find Liquid Concentrations
    C_CO2_v = n_CO2_v / V_v
    C_MEA_v = n_MEA_v / V_v
    C_H2O_v = n_H2O_v / V_v
    C_N2_v = n_N2_v / V_v
    C_O2_v = n_O2_v / V_v
    C_T_v = C_CO2_v + C_MEA_v + C_H2O_v + C_N2_v + C_O2_v

    # Make arrays for DataFrame
    C_v = np.array([C_CO2_v, C_MEA_v, C_H2O_v, C_N2_v, C_O2_v, C_T_v])
    w_v = np.array([w_CO2_v, w_MEA_v, w_H2O_v, w_N2_v, w_O2_v, w_T_v])
    m_v = np.array([m_CO2_v, m_MEA_v, m_H2O_v, m_N2_v, m_O2_v, m_T_v])
    n_v = np.array([n_CO2_v, n_MEA_v, n_H2O_v, n_N2_v, n_O2_v, n_T_v])
    u_v = np.array([u_v, u_v, u_v, u_v, u_v, u_v])
    V_v = np.array([V_v, V_v, V_v, V_v, V_v, V_v])
    y = np.array([y_CO2, y_MEA, y_H2O, y_N2, y_O2, y_T])

    vapor_data = [w_v, m_v, MWs_2, n_v, u_v, V_v, C_v, y]

    vapor_dict = {}
    for i, c in enumerate(columns_v):
        vapor_dict[c] = vapor_data[i]

    vapor = pd.DataFrame(vapor_dict, index=species)
    vapor.style.set_caption("Vapor")

    return vapor, liquid
