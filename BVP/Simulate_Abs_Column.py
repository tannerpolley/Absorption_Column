from scipy.integrate import solve_ivp, odeint
from BVP.ABS_Column import abs_column
from BVP.Solve_BCs import solve_bcs


def simulate_abs_column(inputs, df_param, scales, show_residuals):
    Y_0, shooter_message = solve_bcs(inputs, df_param, scales, show_residuals)

    run_type = 'simulating'

    Fl_z, Fv_0, Tl_z, Tv_0, z, A, P = inputs

    options = {'first_step': 1,
               'max_step': 1,
               'atol': 1e-7,
               }

    options_default = {'max_step': 1,
                       'rtol': .001,
                       'atol': 1e-06,
                       'jac': None,
                       'jac_sparsity': None,
                       'vectorized': False,
                       'first_step': None,

                        'first_step': 1,

                       'atol': 1e-7,
                       }

    # Y = odeint(abs_column, Y_0, z, args=(Fl_z[1], Fv_0[2], Fv_0[3], P, df_param, run_type), tfirst=True).T

    obj = solve_ivp(abs_column,
                    [z[0], z[-1]],
                    Y_0,
                    args=(Fl_z[1], Fv_0[2], Fv_0[3], P, A, df_param, run_type),
                    method='Radau',
                    t_eval=z,
                    vectorized=False,
                    options=options
                    )

    t, Y, nfev, njev, nlu, status, message, success = obj.t, obj.y, obj.nfev, obj.njev, obj.nlu, obj.status, obj.message, obj.success

    return Y, shooter_message
