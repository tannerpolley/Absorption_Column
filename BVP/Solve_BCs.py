from scipy.optimize import root
from BVP.Shooter import shooter


def solve_bcs(inputs, Tl_0_guess, df_param, show_residuals):

    Fl_z, Fv_0, Tl_z, Tv_0, z = inputs[:-1]

    Fv_CO2_0, Fv_H2O_0, Fv_N2_0, Fv_O2_0 = Fv_0

    solve_for_Fl_CO2_0 = False
    solve_for_Fl_H2O_0 = False
    solve_for_Tl_0 = True

    solve_for_list = [solve_for_Fl_CO2_0, solve_for_Fl_H2O_0, solve_for_Tl_0]

    Fl_CO2_0_guess = 3.31034
    Fl_H2O_0_guess = 67.6699

    # with .875 of kv_H2O
    # 323.99596467109725494993
    # 323.99596467109725494992

    # 323.75400802199  Diverges Up
    # 323.75400802198  Collides at 3.2
    # 323.75400802197  Diverges Up
    # 323.75400802196  Collides at 3.24
    # 323.75400802195  Collides at 4.1
    # 323.75400802194  Collides at 4.1
    # 323.75400802193  Doesn't Collide
    # 323.75400802192  Collide at 6? but if H = 9 is used doesn't collide
    # 323.75400802191  Doesn't Collide
    # 323.75400802190  Doesn't Collide
    # 323.75400802180  Doesn't Collide
    # 323.75400802170  Doesn't Collide
    # 323.75400802160  Doesn't Collide
    # 323.75400802100  Doesn't Collide
    # 323.75400802000  Doesn't Collide
    # 323.75400801998  Doesn't Collide
    # 323.75400801997  Doesn't Collide
    # 323.75400801996  Collides at 2.6
    # 323.75400801995  Collides at 2.6
    # 323.75400801990  Collides at 2.6
    # 323.75400801980  Collides at 2.6
    # 323.75400801975  Collides at 2.6
    # 323.75400801974  Collides at 2.6
    # 323.75400801973  Diverges Up
    # 323.75400801972  Diverges Up
    # 323.75400801970  Diverges Up
    # 323.75400801960  Diverges Up
    # 323.75400801950  Diverges Up
    # 323.75400801930  Diverges Up
    # 323.75400801929  Collides at 1.9
    # 323.75400801927  Collides at 1.9
    # 323.75400801925  Collides at 1.9
    # 323.75400801920  Collides at 1.9
    # 323.75400801900  Collides at 1.9
    # 323.75400801000  Collides at 1.9
    # All below this collide earlier and earlier
    # 323.75400805
    # Tl_0_guess = 323.75400802250

    Yv_0_guess = [Fl_CO2_0_guess, Fl_H2O_0_guess, Tl_0_guess]

    # Yv_0_guess = []
    #
    # for i, boolean in enumerate(solve_for_list):
    #     if boolean:
    #         Yv_0_guess.append(Yv_0_poss_guess[i])

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



