from scipy.optimize import root
from BVP.Shooter import shooter
from Parameters import CO2_cap_guess


def solve_bcs(inputs, df_param, show_residuals):

    Fl_0, Fv_z, Tl_0, Tv_z, z = inputs[:-1]

    CO2_cap_guess_dec = CO2_cap_guess*1e-2

    # Top of column Boundary Condition guesses for Vapor Flow Rates and Vapor Temperature
    Fv_CO2_0_guess = Fv_z[0] * (1 - CO2_cap_guess_dec)
    Fv_H2O_0_guess = 3.918873559
    Tv_0_guess = 333

    Yv_0_guess = [Fv_CO2_0_guess, Fv_H2O_0_guess, Tv_0_guess]

    # Prints out the CO2% Captured Guess Value before each run
    CO2_cap_guess_print = (1 - Fv_CO2_0_guess / Fv_z[0]) * 100
    print(f'CO2 Cap Guess/Actual: {CO2_cap_guess_print:.2f}%/')

    shoot = True

    if shoot:

        method = 'Krylov'
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
                           args=(inputs, df_param),
                           method=method,
                           options=options)

        solved_initials, solved, term_msg, n_eval = root_output.x, root_output.success, root_output.message, root_output.nit

        shooter_message = f'Solved? {solved}, with {n_eval:02d} obj function evaluations'

        Fv_CO2_0, Fv_H2O_0, Tv_0 = solved_initials

        Y_0 = [Fl_0[0], Fl_0[2], Fv_CO2_0, Fv_H2O_0, Tl_0, Tv_0]

    else:
        shooter_message = 'No shooting'

        Y_0 = [Fl_0[0], Fl_0[2], Fv_CO2_0_guess, Fv_H2O_0_guess, Tl_0, Tv_0_guess]

    return Y_0, shooter_message



