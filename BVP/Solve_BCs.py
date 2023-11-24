from scipy.optimize import root
from BVP.Shooter import shooter
# from Surrogate import f_CO2_guess, surrogate
from Parameters import Tl_0, Tv_f, CO2_cap_guess
from numpy import append, array
import numpy as np


def solve_bcs(Fl_0, Fv_z, T_0, df_param, temp_dep, show_residuals):

    solve_for_CO2 = True
    solve_for_MEA = False
    solve_for_H2O = True
    solve_for_N2 = False
    solve_for_O2 = False
    solve_for_Tv = False
    solve_for_uv = False

    tests = [solve_for_CO2, solve_for_MEA, solve_for_H2O, solve_for_N2, solve_for_O2, solve_for_Tv, solve_for_uv]

    CO2_cap_guess_dec = CO2_cap_guess*1e-2
    Fv_CO2_0_guess = Fv_z[0] * (1 - CO2_cap_guess_dec)
    Fv_H2O_0_guess = 3.918873559

    Yv_0_guess = [Fv_CO2_0_guess, Fv_H2O_0_guess]
    # Tv_0_guess = 334

    CO2_cap_guess_print = (1 - Fv_CO2_0_guess / Fv_z[0]) * 100
    print(f'CO2 Cap Guess/Actual: {CO2_cap_guess_print:.2f}%/')

    shoot = True

    if shoot:

        method = 'df-sane'
        display = False

        if method == 'df-sane':

            options = {'ftol': 1e-2,
                       'fatol': .25,
                       'maxfev': 50,
                       'line_search': 'cruz',
                       'disp': display,
                       'sigma_0': .1
            }

        elif method == 'Krylov':

            options = {'ftol': .2,
                       'fatol': .25,
                       'maxiter': 50,
                       'disp': display,
                       'line_search': 'armijo',

                       }

        root_output = root(shooter,
                           Yv_0_guess,
                           args=(Fl_0, Fv_z, T_0, df_param, temp_dep),
                           method=method,
                           options=options)

        # bdns = ((.5, 1.5), (0, 100))
        #
        # root_output = minimize(shooter_optimize,
        #                    Yv_0_guess,
        #                    args=(Cl_0, Cv_f, Tl_0, Tv_f, uv_f, Yv_0_known, tests, ul, df_param, temp_dep),
        #                    bounds=bdns)

        solved_initials, solved, term_msg, n_eval = root_output.x, root_output.success, root_output.message, root_output.nit

        shooter_message = f'Solved? {solved}, with {n_eval:02d} obj function evaluations'

        Fv_CO2_0, Fv_H2O_0 = solved_initials

        Y_0 = [Fl_0[0], Fl_0[2], Fv_CO2_0, Fv_H2O_0, T_0[0], T_0[1]]

    else:
        shooter_message = 'No shooting'

        Y_0 = [Fl_0[0], Fl_0[2], Fv_CO2_0_guess, Fv_H2O_0_guess, T_0[0], T_0[1]]

    return Y_0, shooter_message



