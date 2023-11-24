import pandas as pd
from Parameters import MWs, R, Tv_in, P, A, rho_l
import numpy as np


def convert_thesis_data(y_CO2, y_MEA, y_H2O, y_N2, y_O2, m_CO2_v, w_CO2_l, w_MEA_l, w_H2O_l, V_l):

    species = ['CO2', 'MEA', 'H2O', 'N2', 'O2', 'T']
    columns_v = ['w', 'm', 'MW', 'n', 'u', 'C', 'y']
    columns_l = ['w', 'm', 'MW', 'n', 'u', 'C', 'x']

    y_T = sum([y_CO2, y_MEA, y_H2O, y_N2, y_O2])
    ys = [y_CO2, y_MEA, y_H2O, y_N2, y_O2, y_T]

    n = len(species)

    MW = {}
    for i in range(n):
        if i < n - 1:
            MW[species[i]] = MWs[i]
        elif i == n - 1:
            MW[species[i]] = sum(MW[s] for s in species[:-1])

    # Vapor

    # CO2 Conversions
    n_CO2_v = m_CO2_v / MWs[0]
    n_T_v = n_CO2_v / y_CO2
    V_v = n_T_v * R * Tv_in / P  # m3/s

    y = {}
    for i in range(n):

        y[species[i]] = ys[i]

    n_v = {}
    for i in range(n):
        if i < n - 1:
            n_v[species[i]] = n_T_v*y[species[i]]
        elif i == n - 1:
            n_v[species[i]] = n_T_v

    u_v = {}
    for i in range(n):
        u_v[species[i]] = V_v / A

    m_v = {}
    for i in range(n):
        if i < n - 1:
            m_v[species[i]] = round(n_v[species[i]] * MW[species[i]], 5)
        elif i == n - 1:
            m_v[species[i]] = sum(m_v[s] for s in species[:-1])

    C_v = {}
    for i in range(n):
        if i < n - 1:
            C_v[species[i]] = n_v[species[i]] / V_v
        elif i == n - 1:
            C_v[species[i]] = sum(C_v[s] for s in species[:-1])

    w_v = {}
    for i in range(n):
        if i < n - 1:
            w_v[species[i]] = m_v[species[i]] / m_v['T']
        elif i == n - 1:
            w_v[species[i]] = sum(w_v[s] for s in species[:-1])

    vapor_data = [w_v, m_v, MW, n_v, u_v, C_v, y]

    vapor_dict = {}
    for i, c in enumerate(columns_v):
        vapor_dict[c] = vapor_data[i]

    vapor = pd.DataFrame(vapor_dict)

    # Liquid
    m_T_l = V_l * rho_l
    w_ls = np.array([w_CO2_l, w_MEA_l, w_H2O_l, .0, .0, sum([w_CO2_l, w_MEA_l, w_H2O_l])])

    w_l = {}
    for i in range(n):
        w_l[species[i]] = w_ls[i]

    m_l = {}
    for i in range(n):
        if i < n - 1:
            m_l[species[i]] = m_T_l * w_l[species[i]]
        elif i == n - 1:
            m_l[species[i]] = m_T_l

    n_l = {}
    for i in range(n):
        if i < n - 1:
            n_l[species[i]] = m_l[species[i]] / MW[species[i]]
        elif i == n - 1:
            n_l[species[i]] = sum(n_l[s] for s in species[:-1])

    u_l = {}
    for i in range(n):
        u_l[species[i]] = V_l / A

    C_l = {}
    for i in range(n):
        if i < n - 1:
            C_l[species[i]] = n_l[species[i]] / V_l
        elif i == n - 1:
            C_l[species[i]] = sum(C_l[s] for s in species[:-1])

    x = {}
    for i in range(n):
        if i < n - 1:
            x[species[i]] = n_l[species[i]] / n_l['T']
        elif i == n - 1:
            x[species[i]] = sum(x[s] for s in species[:-1])

    # Liquid

    liquid_data = [w_l, m_l, MW, n_l, u_l, C_l, x]

    liquid_dict = {}
    for i, c in enumerate(columns_l):
        liquid_dict[c] = liquid_data[i]

    liquid = pd.DataFrame(liquid_dict)

    return vapor, liquid
