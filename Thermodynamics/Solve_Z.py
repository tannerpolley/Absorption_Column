def solve_Z(A_mix, B_mix, CubicEquationSolver, real):

    m3 = 1
    m2 = B_mix - 1
    m1 = A_mix - 3 * B_mix ** 2 - 2 * B_mix
    m0 = B_mix ** 3 + B_mix ** 2 - A_mix * B_mix

    res = real(CubicEquationSolver.solve(m3, m2, m1, m0))

    try:
        Z_v, Z_l = max(res), min(res)
    except TypeError:
        Z_v, Z_l = .985, .0015

    return Z_v, Z_l
