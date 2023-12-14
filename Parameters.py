import numpy as np
import pandas as pd
from scipy.interpolate import interp1d

# Molecular Weights
MWs_l = np.array([44.01, 61.08, 18.02]) / 1000  # kg/mol
MWs_v = np.array([44.01, 18.02, 28.01, 32]) / 1000  # kg/mol

# Packing Coefficients Mellapak Metal 250Y
a_p = 250  # Packing Surface Area/Volume (m^2/m^3)
Clp = .5  # Packing Coefficient
Chp = .554  # Packing Coefficient
Cvp = .357  # Packing Coefficient
ϵ = .970  # Packing Void Fraction

# Absorption Column Parameters
# 2.008862 Not Clean
# 2.008861 Clean
# with .875 kv_H2O and # 323.99596467109725494993 323.99596467109725494992 as T max and min


H = 6  # Height of Column (m)
D = .64  # Diameter of Column (m)
A = ϵ * np.pi * D ** 2 / 4

n = 101
stages = np.round(np.linspace(0, 1, n))

df_surr = pd.read_csv(r'data\T_surrogate.csv')


def f_interp(name, zi):
    n_int = len(df_surr[name].to_numpy())
    return np.float64(interp1d(np.linspace(0, H, n_int), df_surr[name].to_numpy()[::-1], kind='cubic')(zi))


# Other Constants
g = 9.81  # Gravitational Constant
R = 8.31446  # J/mol-K
# rho_w20 = 966.47 # kg/m3
# mul_w20 = .0010214 # Pa*s
# sigma_w20 = .072967 # N/m

# vapor_species = ['CO2', 'H2O', 'N2', 'O2']


def make_dfs_dict(output_dict, keys_dict, stages):

    sheetnames = list(keys_dict.keys())

    dfs_dict = {}
    for k1 in sheetnames:

        d = {}
        keys = keys_dict[k1]
        array = output_dict[k1]

        for k2, v in zip(keys, array.T):
            d[k2] = v[::-1]

        df = pd.DataFrame(d, index=stages[::-1])
        df.index.name = 'Position'
        # df.index = df.index.astype(int)
        dfs_dict[k1] = df

    return dfs_dict
