from scipy.optimize import root
from BVP.Shooter import shooter


def solve_bcs(inputs, df_param, show_residuals):

    Fl_z, Fv_0, Tl_z, Tv_0, z = inputs[:-1]

    Fv_CO2_0, Fv_H2O_0, Fv_N2_0, Fv_O2_0 = Fv_0

    Fl_CO2_0_guess = 3.31034
    Fl_H2O_0_guess = 67.6699

    # with .875 of kv_H2O
    # 323.99596467109725494993
    # 323.99596467109725494992

    Tl_0_guess = 323.99596467109725494993

    Yv_0_guess = [Fl_CO2_0_guess, Fl_H2O_0_guess, Tl_0_guess]

    shoot = False

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

        Fl_CO2_0, Fl_H2O_0, Tl_0 = solved_initials

        Y_0 = [Fl_CO2_0, Fl_H2O_0, Fv_CO2_0, Fv_H2O_0, Tl_0, Tv_0]

    else:
        shooter_message = 'No shooting'

        Y_0 = [Fl_CO2_0_guess, Fl_H2O_0_guess, Fv_CO2_0, Fv_H2O_0, Tl_0_guess, Tv_0]

    return Y_0, shooter_message



