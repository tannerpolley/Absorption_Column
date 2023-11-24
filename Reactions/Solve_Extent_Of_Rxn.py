
# from numpy import zeros
#
# def solve_extent_of_rxn(Fl, zi):
#
#     ST_1 = [-1, -2, 0, 1, 1, 0]  # 1
#     ST_2 = [-1, -1, -1, 1, 0, 1]  # 2
#
#     Fl_true = zeros(len(ST_1))
#
#     Fl_0 = [Fl[0], Fl[1], Fl[2], 0, 0, 0]
#
#     for i in range(len(Fl_true)):
#         Fl_true[i] = Fl_0[i] + ST_1[i]*f_xi_car(zi) + ST_2[i]*f_xi_bic(zi)
#
#     return Fl_true
