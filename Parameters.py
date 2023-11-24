from numpy import array, pi, linspace, arange, round
import pandas as pd
from scipy.interpolate import interp1d

print_progress = False

if print_progress:

    print_progress_temp = True
    print_progress_flux = True
    print_progress_thermo = True
    print_progress_heat = True

else:

    print_progress_temp = False
    print_progress_flux = False
    print_progress_thermo = False
    print_progress_heat = False


# Molecular Weights
MWs_l = array([44.01, 61.08, 18.02]) / 1000 # kg/mol
MWs_v = array([44.01, 18.02, 28.01, 32]) / 1000 # kg/mol

T_AB = 319.15  # Temperature of the system (K)
Tl_0 = 315.39    # Temperature of Liquid into Abs (K)
Tv_f = 319.15    # Temperature of Vapor into Abs (K)
T_amb = Tv_f
Tv_in = Tv_f
CO2_cap_guess = 90.4

# Packing Coefficients Mellapak 252Y
# a_T = 250  # Packing Surface Area/Volume (m^2/m^3)
# Clp = .20333  # Packing Coefficient
# Chp = .544  # Packing Coefficient
# Cvp = .352251  # Packing Coefficient
# eps = .980  # Packing Void Fraction

# Packing Coefficients Mellapak Metal 250Y
a_p = 250  # Packing Surface Area/Volume (m^2/m^3)
Clp = .5  # Packing Coefficient
Chp = .554  # Packing Coefficient
Cvp = .357  # Packing Coefficient
ϵ = .970  # Packing Void Fraction

# a_p = 200  # Packing Surface Area/Volume (m^2/m^3)
# Clp = .971  # Packing Coefficient
# Chp = .547  # Packing Coefficient
# Cvp = .390  # Packing Coefficient
# ϵ = .979  # Packing Void Fraction

# Absorption Column Parameters
P = 109180  # (Pa) Pressure of the System
H = 6  # Height of Column (m)
D = .64  # Diameter of Column (m)
A = ϵ * pi * D ** 2 / 4

n = 91
# dz = .2
z = linspace(0, H, n)

df = pd.read_csv(r'C:\Users\364265\Documents\Absorption_Column_7.8 - Shooting Method with Energy dependence but with fudge factors\data\T_surrogate.csv')

f_Tl = interp1d(z,  df['Tl'].to_numpy(), kind='cubic')
f_Tv = interp1d(z,  df['Tv'].to_numpy(), kind='cubic')
f_E = interp1d(z, df['E'].to_numpy(), kind='cubic')

f_duvdz = interp1d(z, df['duvdz'].to_numpy(), kind='cubic')
# f_kl_CO2 = interp1d(z, df['kl_CO2'].to_numpy(), kind='cubic')
# f_kv_CO2 = interp1d(z, df['kv_CO2'].to_numpy(), kind='cubic')
# f_kv_H2O = interp1d(z, df['kv_H2O'].to_numpy(), kind='cubic')
f_DF_H2O = interp1d(z, df['DF_H2O'].to_numpy(), kind='cubic')
# f_Dv_CO2 = interp1d(z, df['Dv_CO2'].to_numpy(), kind='cubic')
# f_Dv_H2O = interp1d(z, df['Dv_H2O'].to_numpy(), kind='cubic')
# f_muv_mix = interp1d(z, df['mu_v'].to_numpy(), kind='cubic')
f_Hv_CO2 = interp1d(z, df['Hv_CO2'].to_numpy(), kind='cubic')
f_Hv_H2O = interp1d(z, df['Hv_H2O'].to_numpy(), kind='cubic')


f_Cl_CO2 = interp1d(z, df['Cl_CO2'].to_numpy(), kind='cubic')
f_Cl_MEA = interp1d(z, df['Cl_MEA'].to_numpy(), kind='cubic')
f_Cl_H2O = interp1d(z, df['Cl_H2O'].to_numpy(), kind='cubic')
f_Cl_MEAH = interp1d(z, df['Cl_MEAH'].to_numpy(), kind='cubic')
f_Cl_MEACOO = interp1d(z, df['Cl_MEACOO'].to_numpy(), kind='cubic')
f_Cl_HCO3 = interp1d(z, df['Cl_HCO3'].to_numpy(), kind='cubic')

f_UT = interp1d(z, df['UT'].to_numpy(), kind='cubic')

def f_interp(name, zi):
    return interp1d(z, df[name].to_numpy(), kind='cubic')(zi)





# z = arange(0, H+dz, dz)
# n = len(z)
stages = linspace(0, n, n)
stages = round(arange(0, n, 1), 0)

# dz = round(z[1] - z[0], 3)

# Other Constants
g = 9.81  # Gravitational Constant
R = 8.31446  # J/mol-K
# rho_l = 1000  # kg/m3 Liquid Inlet Density (Like Water)
rho_l = 1022.4252
rho_w20 = 966.47 # kg/m3
mul_w20 = .0010214 # Pa*s
sigma_w20 = .072967 # N/m


def f_w_MEA(z, α, x_MEA):
    return x_MEA - (1 + α + (MWs_l[1] / MWs_l[2]) * ((1 - z) / z)) ** (-1)


# Interaction Parameters
CO2_MEA = .16
CO2_H2O = .065
CO2_N2 = -.0149
CO2_O2 = -.04838
MEA_H2O = -.18
MEA_N2 = 0
MEA_O2 = 0
H2O_N2 = 0
H2O_O2 = 0
N2_O2 = -.00978

# k_int = array([[0,       CO2_MEA, CO2_H2O, CO2_N2, CO2_O2],
#                [CO2_MEA,       0, MEA_H2O, MEA_N2, MEA_O2],
#                [CO2_H2O, MEA_H2O,       0, H2O_N2, H2O_O2],
#                [CO2_N2,   MEA_N2,  H2O_N2,      0, N2_O2],
#                [CO2_O2,   MEA_O2,  H2O_O2,  N2_O2,     0]])

k_int = array([[0,       CO2_MEA, CO2_H2O],
               [CO2_MEA,       0, MEA_H2O],
               [CO2_H2O, MEA_H2O,       0]])

# EOS Parameters
data_EOS = array([[304.21, 7.383*1e6, 0.228],
                  [678.2, 4.45*1e6, 0.864],
                  [647.096, 22.12*1e6, 0.344],
                  [126.2, 3.394*1e6, 0.040],
                  [154.58, 5.043*1e6, 0.022]])

data_EOS = data_EOS.T

data = data_EOS

Tc, Pc, ac = data

liquid_species = ['CO2', 'MEA', 'H2O']
vapor_species = ['CO2', 'H2O', 'N2', 'O2']


Fl_keys = ['Fl_CO2', 'Fl_MEA', 'Fl_H2O']
Fl_true_keys = ['Fl_CO2_i', 'Fl_MEA_i', 'Fl_H2O_i', 'Fl_MEAH_i', 'Fl_MEACOO_i', 'Fl_HCO3_i']
Fv_keys = ['Fv_CO2', 'Fv_H2O', 'Fv_N2', 'Fv_O2']
Cl_keys = ['Cl_CO2', 'Cl_MEA', 'Cl_H2O']
Cl_true_keys = ['Cl_CO2_i', 'Cl_MEA_i', 'Cl_H2O_i', 'Cl_MEAH_i', 'Cl_MEACOO_i', 'Cl_HCO3_i']
Cv_keys = ['Cv_CO2', 'Cv_MEA', 'Cv_H2O', 'Cv_N2', 'Cv_O2']

y_keys = ['y_CO2', 'y_H2O', 'y_N2', 'y_O2']
x_keys = ['x_CO2', 'x_MEA', 'x_H2O']
x_true_keys = ['x_CO2_i', 'x_MEA_i', 'x_H2O_i', 'x_MEAH_i', 'x_MEACOO_i', 'x_HCO3_i']
flux_keys = ['N_CO2', 'N_H2O']
ql_keys = ['ql_CO2', 'ql_H2O', 'q_trn', 'q_l', 'dTl_dz']
qv_keys = ['qv_CO2', 'qv_H2O', 'q_trn', 'q_v', 'dTv_dz']
hydr_keys = ['ul', 'uv', 'h_L', 'a_e', 'a_e*A']

PEQ_keys =  ['DF_CO2', 'DF_H2O', 'Pv_CO2', 'Pl_CO2', 'Pv_H2O', 'Pl_H2O', 'H_CO2_mix', 'Psat_H2O']
k_mx_keys = ['kl_CO2', 'kv_CO2', 'kv_H2O', 'KH', 'E', 'k_El', 'k_Ev']
T_keys = ['T_l', 'T_l2', 'T_v', 'T_v2']
liquid_prop_keys = ['rho_mol', 'rho_mass', 'mu', 'sigma', 'Dl_CO2']
vapor_prop_keys = ['rho_mol', 'rho_mass', 'mu', 'Dv_CO2', 'Dv_H2O']
Cpl_keys = ['Cpl_CO2', 'Cpl_MEA', 'Cpl_H2O']
Cpv_keys = ['Cpv_CO2', 'Cpv_H2O', 'Cpv_N2', 'Cpv_O2']

keys_dict = {'Fl': Fl_keys + Fl_true_keys,
             'Fv': Fv_keys,
             'Cl': Cl_keys + Cl_true_keys,
             'Cv': Cv_keys,
             'x': x_keys + x_true_keys,
             'y': y_keys,
             'T': T_keys,
             'PEQ': PEQ_keys,
             'hydr': hydr_keys,
             'k_mx': k_mx_keys,
             'N': flux_keys,
             'Prop_l': liquid_prop_keys + Cpl_keys,
             'Prop_v': vapor_prop_keys + Cpv_keys,
             'ql': ql_keys,
             'qv': qv_keys}

sheetnames = list(keys_dict.keys())


def make_dfs_dict(output_dict):

    dfs_dict = {}
    for k1 in sheetnames:

        d = {}
        keys = keys_dict[k1]
        array = output_dict[k1]

        for k2, v in zip(keys, array.T):
            d[k2] = v

        df = pd.DataFrame(d, index=stages)
        df.index.name = 'Stages'
        df.index = df.index.astype(int)
        dfs_dict[k1] = df

    return dfs_dict
