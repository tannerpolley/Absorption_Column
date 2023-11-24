def solve_Phi(Tl, Tv, x, args1, args2):

    Tc, Pc, ac, P, k_int, R, solve_Z, CubicEquationSolver, print_progress_thermo, print_progress_temp, flag = args1
    log, column_stack, sum, real, array, exp, sqrt, append = args2

    Tr = Tl / Tc[:3]

    κ = 0.37464 + 1.54226 * ac[:3] - 0.26992 * ac[:3] ** 2

    α = (1 + κ * (1 - Tr ** .5)) ** 2

    a = α * (.45724 * R ** 2 * Tc[:3] ** 2 / Pc[:3])
    aj = a.T
    b = .0778 * R * Tc[:3] / Pc[:3]

    col_a = sqrt(a * aj[0]) * (1 - k_int[0])
    col_b = sqrt(a * aj[1]) * (1 - k_int[1])
    col_c = sqrt(a * aj[2]) * (1 - k_int[2])

    a_ij = column_stack([col_a, col_b, col_c])

    col_a = a_ij[0] * x * x[0]
    col_b = a_ij[1] * x * x[1]
    col_c = a_ij[2] * x * x[2]

    xa_ij = column_stack([col_a, col_b, col_c])

    a_mix = sum(xa_ij)
    b_mix = sum(b * x)

    A_mix = a_mix * P / (R ** 2 * Tl ** 2)
    B_mix = b_mix * P / (R * Tl)

    col_1 = sum(a_ij[:, 0] * x)
    col_2 = sum(a_ij[:, 1] * x)
    col_3 = sum(a_ij[:, 2] * x)

    σ = array([col_1, col_2, col_3])

    Z_v, Z_l = solve_Z(A_mix, B_mix, CubicEquationSolver, real)

    def f_Φ(phase):
        if phase == 'liquid':
            Z = Z_l
        elif phase == 'vapor':
            Z = Z_v

        # inside = (b / b_mix * (Z - 1) - log((Z - B_mix)) -
        #       ((A_mix / (2.828427 * B_mix) * (2 * σ / a_mix - b / b_mix)) *
        #       log(((Z + 2.414 * B_mix) / (Z - .414 * B_mix)))))

        term_1 = b / b_mix * (Z - 1)
        # print(b_mix, P, R, Tl, B_mix)
        if Z - B_mix < 0:
            # print(b_mix, P, R, Tl, B_mix, 'error')
            # print(B_mix, x)
            # print((Z - B_mix))
            # print(Z - .414 * B_mix)
            term_2 = log(abs(Z - B_mix))
        else:
            term_2 = log((Z - B_mix))
        term_3 = A_mix / (2.828427 * B_mix)
        term_4 = (2 * σ / a_mix - b / b_mix)
        term_5 = (Z + 2.414 * B_mix)
        term_6 = (Z - .414 * B_mix)

        # if term_5 < 0:
        #     print(x, 'Error')
        #     return 'Error'

        term_7 = log(term_5 / term_6)
        term_8 = (term_3 * term_4 * term_7)

        inside = term_1 - term_2 - term_8

        Φ = exp(inside)

        return Φ

    Φ_v = f_Φ('vapor')
    Φ_l = f_Φ('liquid')

    Φ_v = append(Φ_v, [1, 1])
    Φ_l = append(Φ_l, [1, 1])

    Z_arr = array([Z_v, Z_l])

    return Z_arr, [Φ_v, Φ_l]
