from scipy.integrate import solve_ivp
from BVP.ABS_Column import abs_column


def shooter(X, inputs, df_param):

    Fl_CO2_0, Fl_H2O_0, Tl_0 = X

    Fl_z, Fv_0, Tl_z, Tv_0, z, P = inputs

    Fl_CO2_z, Fl_MEA_z, Fl_H2O_z = Fl_z
    Fv_CO2_0, Fv_H2O_0, Fv_N2_0, Fv_O2_0 = Fv_0

    Y_0 = [Fl_CO2_0, Fl_H2O_0, Fv_CO2_0, Fv_H2O_0, Tl_0, Tv_0]

    run_type = 'shooting'

    Y = solve_ivp(abs_column,
                  [z[0], z[-1]],
                  Y_0,
                  args=(Fl_MEA_z, Fv_N2_0, Fv_O2_0, P, df_param, run_type),
                  method='Radau').y

    Fl_CO2_z_sim, Fl_H2O_z_sim, Tl_z_sim = Y[0, -2], Y[1, -2], Y[4, -2]

    eq1 = Fl_CO2_z_sim - Fl_CO2_z
    eq2 = Fl_H2O_z_sim - Fl_H2O_z
    eq3 = Tl_z_sim - Tl_z

    eqs = [eq1, eq2, eq3]

    return eqs
