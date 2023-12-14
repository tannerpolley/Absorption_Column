from scipy.integrate import solve_ivp, odeint
from BVP.ABS_Column import abs_column
from BVP.Solve_BCs import solve_bcs


def simulate_abs_column(inputs, df_param, show_residuals):

    Y_0, shooter_message = solve_bcs(inputs, df_param, show_residuals)

    if shooter_message == 'This will break the model':
        return -1000, shooter_message

    run_type = 'simulating'

    Fl_z, Fv_0, Tl_0, Tv_z, z, P = inputs

    obj = solve_ivp(abs_column,
                    [z[0], z[-1]],
                    Y_0,
                    args=(Fl_z[1], Fv_0[2], Fv_0[3], P, df_param, run_type),
                    method='Radau',
                    t_eval=z,
                    )

    t, Y, nfev, njev, nlu, status, message, success = obj.t, obj.y, obj.nfev, obj.njev, obj.nlu, obj.status, obj.message, obj.success

    return Y, shooter_message
