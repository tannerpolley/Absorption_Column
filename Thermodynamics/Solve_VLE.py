from numpy import log, column_stack, sum, real, array, exp, sqrt, append, linalg
from scipy.optimize import root
from Thermodynamics.Solve_Interface_xy import solve_interface_xy
from Properties.Henrys_Law import henrys_law
from Thermodynamics.Solve_Phi import solve_Phi
from Thermodynamics.Solve_Z import solve_Z
from Parameters import Tc, Pc, ac, P, k_int, R, print_progress_thermo, print_progress_temp
from misc import CubicEquationSolver


def solve_vle(Cl, Cv, Tl, Tv, df_param):
    Cl_T = sum(Cl)
    Cv_T = sum(Cv)

    y = Cv / Cv_T
    x = Cl / Cl_T

    y_CO2, y_MEA, y_H2O, y_N2, y_O2 = y
    x_CO2, x_MEA, x_H2O, x_N2, x_O2 = x[:5]

    guesses = array([x_CO2, y_MEA, y_H2O])

    args1 = [Tc, Pc, ac, P, k_int, R, solve_Z, CubicEquationSolver, print_progress_thermo,
             print_progress_temp, 'in']

    args2 = [log, column_stack, sum, real, array, exp, sqrt, append]

    try:
        xi_CO2, yi_MEA, yi_H2O = root(solve_interface_xy,
                                      guesses,
                                      method='df-sane',
                                      args=(Tl, Tv, x, y, solve_Phi, args1, args2)).x
    except ValueError:
        print('hi')
        print(Cl)
        print(Cv)
        solve_vle(Cl, Cv, Tl, Tv, df_param)

    xi_N2, xi_O2 = 0, 0

    xi = xi_CO2, x_MEA, x_H2O

    flag = 'out'
    Z, Φ = solve_Phi(Tl, Tv, xi, args1, args2)
    K = Φ[1] / Φ[0]

    Cl_CO2, Cl_MEA, Cl_H2O, Cl_N2, Cl_O2 = Cl[:5]
    Cv_CO2, Cv_MEA, Cv_H2O, Cv_N2, Cv_O2 = Cv

    A = array([[1 - xi_CO2, 0, 0, -xi_CO2, -xi_CO2],
               [0, 1 - yi_MEA, -yi_MEA, 0, 0],
               [0, -yi_H2O, 1 - yi_H2O, 0, 0],
               [-xi_N2, 0, 0, 1 - xi_N2, -xi_N2],
               [-xi_O2, 0, 0, -xi_O2, 1 - xi_CO2]])

    B = array([xi_CO2 * (Cl_MEA + Cl_H2O),
               yi_MEA * (Cv_CO2 + Cv_O2 + Cv_N2),
               yi_H2O * (Cv_CO2 + Cv_O2 + Cv_N2),
               xi_N2 * (Cl_MEA + Cl_H2O),
               xi_O2 * (Cl_MEA + Cl_H2O)])

    Cli_CO2, Cvi_MEA, Cvi_H2O, Cli_N2, Cli_O2 = linalg.solve(A, B)

    Cli = array([Cli_CO2, Cl_MEA, Cl_H2O, Cli_N2, Cli_O2])
    Cvi = array([Cv_CO2, Cvi_MEA, Cvi_H2O, Cv_N2,  Cv_O2])

    return Cli, Cvi, K, Z, Φ
