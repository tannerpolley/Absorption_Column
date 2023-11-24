from numpy import array, sum
from scipy.optimize import root
from Transport.Solve_MassTransfer import solve_masstransfer
from Thermodynamics.Solve_Driving_Force import solve_driving_force
from Transport.Solve_Flux import solve_flux
from Reactions.Solve_ChemEQ import solve_ChemEQ_log
from Transport.Solve_Temperature_Change import solve_temperature_change
from Properties.Henrys_Law import henrys_law
from Properties.Density import liquid_density, vapor_density
from Properties.Viscosity import viscosity
from Properties.Surface_Tension import surface_tension
from Properties.Diffusivity import liquid_diffusivity, vapor_diffusivity
from Properties.Heat_Capacity import heat_capacity
from Properties.Thermal_Conductivity import thermal_conductivity
from Parameters import f_w_MEA, f_Tl, f_Tv, A, P, MWs_l
from Parameters import f_Cl_CO2, f_Cl_MEA, f_Cl_H2O,f_Cl_MEAH,f_Cl_MEACOO, f_Cl_HCO3, \
                        f_UT, f_Hv_CO2, f_Hv_H2O, f_duvdz, f_interp
from Transport.Heat_Transfer import heat_transfer


def abs_column(zi, Y, Fl_MEA, Fv_N2, Fv_O2, df_param, temp_dep, run_type):

    Fl_CO2, Fl_H2O, Fv_CO2, Fv_H2O, Tl_2, Tv_2 = Y

    Tl = f_Tl(zi)
    Tv = f_Tv(zi)

    Tl = Tl_2
    Tv = Tv_2

    Fl_T = Fl_CO2 + Fl_MEA + Fl_H2O
    Fv_T = Fv_CO2 + Fv_H2O + Fv_N2 + Fv_O2

    Fl = [Fl_CO2, Fl_MEA, Fl_H2O]
    Fv = [Fv_CO2, Fv_H2O, Fv_N2, Fv_O2]

    x = [Fl[i]/Fl_T for i in range(len(Fl))]
    y = [Fv[i] / Fv_T for i in range(len(Fv))]

    w = [MWs_l[i]*x[i] / sum([MWs_l[j]*x[j] for j in range(len(Fl))]) for i in range(len(Fl))]

    alpha = x[0] / x[1]

    w_MEA = root(f_w_MEA, array([.3]), args=(alpha, x[1])).x[0]

    # -------------------------------- Properties ---------------------------------------------

    # Density
    rho_mol_l, rho_mass_l = liquid_density(Tl, w_MEA, alpha, df_param)
    rho_mol_v, rho_mass_v = vapor_density(Tv, P, y)

    # Viscosity
    muv_mix, mul_mix, mul_H2O = viscosity(Tl, Tv, y, w_MEA, alpha, df_param)

    # Surface Tension
    sigma = surface_tension(Tl, x, w_MEA, alpha, df_param)

    # Diffusion
    Dl_CO2 = liquid_diffusivity(Tl, rho_mol_l*x[1], df_param)
    Dv_CO2, Dv_H2O, Dv_N2, Dv_O2, Dv_T = vapor_diffusivity(Tv, y, P, df_param)

    # Heat Capacity

    Cpl = heat_capacity(Tl, 'liquid', x, w)
    Cpv = heat_capacity(Tv, 'vapor', x, w)

    Cpv_T = sum([Cpv[i]*y[i] for i in range(len(y))])

    kt_CO2 = thermal_conductivity(Tv, 'CO2', 'vapor')

    H_CO2_mix = henrys_law(x, Tl)

    # ------------------------------ Chemical Equilibrium --------------------------------------

    Cl = [x[i] * rho_mol_l for i in range(len(x))]
    Cv = [y[i] * rho_mol_v for i in range(len(x))]

    guesses = f_Cl_CO2(zi), f_Cl_MEA(zi), f_Cl_H2O(zi), \
        f_Cl_MEAH(zi), f_Cl_MEACOO(zi), f_Cl_HCO3(zi)

    Cl_true = solve_ChemEQ_log(Cl, Tl, guesses)

    x_true = [Cl_true[i] / rho_mol_l for i in range(len(Cl_true))]

    Fl_true = [x_true[i] * Fl_T for i in range(len(x_true))]

    # ------------------------------ Transport --------------------------------------

    ul = Fl_T / (A * rho_mol_l)
    uv = Fv_T / (A * rho_mol_v)

    kv_CO2, kv_H2O, kv_T, KH, k_mxs, uv, Dv_T, a_e, hydr = solve_masstransfer(rho_mass_l, rho_mass_v, mul_mix, muv_mix,
                                                                              sigma, Dl_CO2, Dv_CO2, Dv_H2O, Dv_T, Tv,
                                                                              ul, uv, H_CO2_mix, zi)
    UT = heat_transfer(P, kv_CO2, kt_CO2, Cpv_T, rho_mol_v, Dv_CO2)

    # ------------------------------ Thermodynamics --------------------------------------

    DF_CO2, DF_H2O, PEQ = solve_driving_force(Cl_true, x_true, y, Tl, Tv, KH, H_CO2_mix, zi)

    N_CO2 = kv_CO2 * DF_CO2
    N_H2O = kv_H2O * DF_H2O

    # Mass Balance

    dFl_CO2_dz = N_CO2*a_e*A
    dFl_H2O_dz = N_H2O*a_e*A

    dFv_CO2_dz = N_CO2*a_e*A
    dFv_H2O_dz = N_H2O*a_e*A

    Sigma_Fl_Cpl = sum([Fl[i] * Cpl[i] for i in range(len(Fl))])
    Sigma_Fv_Cpv = sum([Fv[i] * Cpv[i] for i in range(len(Fv))])

    Hl_CO2 = 84000
    Hl_H2O = 44000

    ql_CO2 = N_CO2*Hl_CO2
    ql_H2O = N_H2O*Hl_H2O
    # q_trn = f_UT(zi) * (Tv - (Tl + f_interp('Tl_fudge', zi)*.8))
    q_trn = UT * (Tv - Tl)
    # q_trn_2 = f_UT(zi) * (Tv - (Tl + f_interp('Tl_fudge', zi)))
    # q_trn = f_interp('q_trn', zi)

    qv_CO2 = N_CO2*f_Hv_CO2(zi)
    qv_H2O = N_H2O*f_Hv_H2O(zi)

    ql =  ql_CO2 + ql_H2O + q_trn
    qv = q_trn

    k_El = a_e*A/Sigma_Fl_Cpl
    k_Ev = a_e*A/Sigma_Fv_Cpv

    dTl_dz = ql*a_e*A/Sigma_Fl_Cpl
    dTv_dz = qv*a_e*A/Sigma_Fv_Cpv

    diffeqs = [dFl_CO2_dz, dFl_H2O_dz, dFv_CO2_dz, dFv_H2O_dz, dTl_dz, dTv_dz]

    if run_type == 'saving':

        output_dict = {'Fl': Fl + list(Fl_true),
                       'Fv': Fv,
                       'Cl': Cl + list(Cl_true),
                       'Cv': Cv,
                       'x': x + list(x_true),
                       'y': y,
                       'T': [Tl, Tl_2, Tv, Tv_2],
                       'PEQ': PEQ,
                       'hydr': hydr + [a_e*A],
                       'k_mx': k_mxs + [k_El, k_Ev],
                       'N': [N_CO2, N_H2O],
                       'Prop_l': [rho_mol_l, rho_mass_l, mul_mix, sigma, Dl_CO2] + Cpl,
                       'Prop_v': [rho_mol_v, rho_mass_v, muv_mix, Dv_CO2, Dv_H2O] + Cpv,
                       'ql': [ql_CO2, ql_H2O, q_trn, ql, dTl_dz],
                       'qv': [qv_CO2, qv_H2O, q_trn, qv, dTv_dz],
                       }
        return output_dict

    return diffeqs


