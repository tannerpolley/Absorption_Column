from scipy.integrate import solve_ivp
from BVP.ABS_Column import abs_column


def shooter(X, inputs, df_param):

    Fv_CO2_0, Fv_H2O_0, Tv_0 = X

    Fl_0, Fv_z, Tl_0, Tv_z, z, P = inputs

    Fl_CO2_0, Fl_MEA_0, Fl_H2O_0 = Fl_0
    Fv_CO2_z, Fv_H2O_z, Fv_N2_z, Fv_O2_z = Fv_z

    Y_0 = [Fl_CO2_0, Fl_H2O_0, Fv_CO2_0, Fv_H2O_0, Tl_0, Tv_0]

    run_type = 'shooting'

    Y = solve_ivp(abs_column,
                  [z[0], z[-1]],
                  Y_0,
                  args=(Fl_MEA_0, Fv_N2_z, Fv_O2_z, P, df_param, run_type),
                  method='Radau').y

    Fv_CO2_z_sim, Fv_H2O_z_sim, Tv_z_sim = Y[2, -1], Y[3, -1], Y[-1, -1]

    eq1 = Fv_CO2_z_sim - Fv_CO2_z
    eq2 = Fv_H2O_z_sim - Fv_H2O_z
    eq3 = Tv_z_sim - Tv_z

    eqs = [eq1, eq2, eq3]

    return eqs
