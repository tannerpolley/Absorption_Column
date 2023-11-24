from scipy.integrate import solve_ivp
from BVP.ABS_Column import abs_column
from Parameters import z
from BVP.Solve_BCs import solve_bcs


def simulate_abs_column(Fl_0, Fv_z, T_0, df_param, temp_dep, show_residuals):

    n_l = len(Fl_0)
    Y_0, shooter_message = solve_bcs(Fl_0, Fv_z, T_0,
                                     df_param,
                                     temp_dep,
                                     show_residuals)

    if shooter_message == 'This will break the model':
        return -1000, shooter_message

    run_type = 'simulating'


    obj = solve_ivp(abs_column,
                    [z[0], z[-1]],
                    Y_0,
                    args=(Fl_0[1], Fv_z[2], Fv_z[3], df_param, temp_dep, run_type),
                    method='Radau',
                    t_eval=z
                    )

    t, Y, nfev, njev, nlu, status, message, success = obj.t, obj.y, obj.nfev, obj.njev, obj.nlu, obj.status, obj.message, obj.success

    return Y, shooter_message
