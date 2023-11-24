def solve_interface_xy(z, Tl, Tv, x, y, solve_Phi, args1, args2):

    y_CO2, y_MEA, y_H2O, y_N2, y_O2 = y
    x_CO2, x_MEA, x_H2O, x_N2, x_O2 = x[:5]

    xi_CO2, yi_MEA, yi_H2O = z
    xi = xi_CO2, x_MEA, x_H2O

    flag = 'in'

    try:
        Z, Φ = solve_Phi(Tl, Tv, xi, args1, args2)
    except:
        return 'Error'

    K = Φ[1] / Φ[0]

    K_CO2, K_MEA, K_H2O, K_N2, K_O2 = K
    eq1 = xi_CO2 - y_CO2 / K_CO2
    eq2 = yi_MEA - K_MEA * x_MEA
    eq3 = yi_H2O - K_H2O * x_H2O
    eqs = [eq1, eq2, eq3]

    return eqs

