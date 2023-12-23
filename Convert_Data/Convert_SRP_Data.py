from Parameters import MWs_l
import numpy as np


def convert_SRP_data(X, n):

    Tl_z, Tv_0, Fl_z, Fv_0, alpha, y_H2O, y_CO2, H = X

    w_MEA = .325
    P = 107560
    D = .43
    A = np.pi * D ** 2 / 4
    z = np.linspace(0, H, n)

    alpha_O2_N2 = 0.08485753604

    MW_MEA = MWs_l[1]
    MW_H2O = MWs_l[2]

    x_MEA = ((1 + alpha + (MW_MEA/MW_H2O))*(1-w_MEA)/w_MEA)**-1
    x_CO2 = x_MEA*alpha
    x_H2O = 1 - x_CO2 - x_MEA

    Fl_CO2_z = x_CO2*Fl_z
    Fl_MEA_z = x_MEA*Fl_z
    Fl_H2O_z = x_H2O*Fl_z

    y_N2 = (1 - y_CO2 - y_H2O)/(1 + alpha_O2_N2)
    y_O2 = alpha_O2_N2*y_N2

    Fv_CO2_0 = Fv_0 * y_CO2
    Fv_H2O_0 = Fv_0 * y_H2O
    Fv_N2_0 = Fv_0 * y_N2
    Fv_O2_0 = Fv_0 * y_O2

    Fl = Fl_CO2_z, Fl_MEA_z, Fl_H2O_z
    Fv = Fv_CO2_0, Fv_H2O_0, Fv_N2_0, Fv_O2_0
    return Fl, Fv, Tl_z, Tv_0, z, A, P



