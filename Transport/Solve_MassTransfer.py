from Parameters import a_p, Clp, Chp, Cvp, ϵ, g, R, f_interp, A
import numpy as np
log = np.log
exp = np.exp


def solve_masstransfer(rho_mass_l, rho_mass_v, mul_mix, muv_mix, sigma, Dl_CO2, Dv_CO2, Dv_H2O, Dv_T, Tl, Tv, ul, uv, H_CO2_mix, Cl_MEA_true, zi):

    # Compute wetted hydraulic specific transfer area
    Re = (ul * rho_mass_l) / (a_p * mul_mix)
    if Re < 5:
        a_h = a_p * Chp * (rho_mass_l * ul / (a_p * mul_mix)) ** .15 * (a_p * ul ** 2 / 9.81) ** .1
    else:
        a_h = a_p * .85 * Chp * (rho_mass_l * ul / (a_p * mul_mix)) ** .25 * (a_p * ul ** 2 / 9.81) ** .1


    # Liquid Hold Up at the loading point
    h_L_s = (12 * ul * a_p ** 2 * mul_mix / (9.81 * rho_mass_l)) ** (1 / 3) * (a_h / a_p) ** (2 / 3)
    #
    # # Liquid Hold Up at the flooding point
    # h_L_Fl = 2.2*h_L_s*(mul_mix*rho_w20/rho_l/mul_w20)**(.05)
    #
    # F_LG = L/G*(rho_v/rho_l)
    #
    # if F_LG <= .4:
    #     n_Fl = -.194
    #     C_Fl = 2.339
    # elif F_LG > .4:
    #     n_Fl = -.708
    #     C_Fl = 2.339*.6244*(mul_mix/muv_mix)**(.1028)
    #
    # # Resistance Coefficient
    # Psi_Fl = g/C_Fl**2*(F_LG*(mul_mix/muv_mix)**.2)**(-2*n_Fl)
    #
    # # Flooding Velocity
    # uv_Fl = (2 * g / Psi_Fl) ** .5 * ((ϵ - h_L_Fl) ** 1.5 / ϵ ** .5) * (h_L_Fl / a_p) ** .5 * (rho_l / rho_v) ** .5
    #
    # Fl_frc = uv/uv_Fl
    #
    # h_L = h_L_s + (h_L_Fl - h_L_s)*(Fl_frc)**13

    Apar = 11.4474
    B = 0.6471
    Alpha = 3.185966

    h_L = np.exp(log(Apar) + B * (log(Alpha) + log(ul) + 1/3*(log(mul_mix) - log(rho_mass_l))))
    h_V = ϵ - h_L
    # h_L = h_L_s

    uv = uv / (1 - h_L)

    if h_L > ϵ:
        print('Error: Overall Vapor Mass Transfer Coefficient will get an error')

    d_h = 4 * ϵ / a_p
    Lp = A*a_p/ϵ

    # Compute interfacial area
    ## Interfacial Area of loading point
    ae_ap_s = 1.5 * (a_p/d_h) ** (-.5) * (rho_mass_l * ul * d_h / mul_mix) ** (-.2) * (rho_mass_l * ul * d_h / sigma) ** .75 * \
         (ul ** 2 / (g * d_h)) ** (-.45)
    #
    # ## Interfacial Area of Flooding point
    # ae_ap_Fl = 7*(sigma/sigma_w20)**.56 * ae_ap_s
    #
    # ## Interfacial Area
    # ae_ap = ae_ap_s + (ae_ap_Fl - ae_ap_s) * (uv / uv_Fl) ** 13
    # a_e = ae_ap*a_p

    a_e = np.exp(log(a_p) + log(1.42) + .116*(log(rho_mass_l) - log(sigma) + (1 / 3) * log(g) + (4 / 3) * (log(ul) + log(A) - log(Lp))))
    # a_e = ae_ap_s

    def f_kl(Dl):
        kl = exp(log(Clp) + (1/6) * log(12) + .5 * (log(ul) + log(Dl) - log(h_L) - log(d_h)))
        # kl = Clp * (rho_mass_l * g / mul_mix) ** (1 / 6) * (a_p * Dl / (4 * ϵ)) ** (1 / 2) * (ul / a_p) ** (1 / 3)
        return kl

    def f_kv(Dv):

        kv = exp(log(Cvp) - log(R) - log(Tv) \
             - 0.5 * log(h_V) \
             + 0.5 * (log(a_p) - log(d_h)) \
             + (2 / 3) * log(Dv) \
             + (1 / 3) * (log(muv_mix) - log(rho_mass_v)) \
             + (3 / 4) * (log(uv) + log(rho_mass_v) - log(a_p) - log(muv_mix)))
        # print(log(Cvp), log(R), log(Tv), log(h_V), log(a_p), log(d_h), log(Dv), log(muv_mix), log(rho_mass_v), log(uv))
        return kv

    kv_CO2, kv_H2O = f_kv(Dv_CO2), f_kv(Dv_H2O)
    kv_T = f_kv(Dv_T)*(R*Tv)

    kl_CO2 = f_kl(Dl_CO2)

    E = f_interp('E', zi)
    # kl_CO2 = f_interp('kl_CO2', zi)
    # kv_CO2 = f_interp('kv_CO2', zi)
    # kv_H2O = f_interp('kv_H2O', zi)
    k2 = 3.1732e9*exp(-4936.6/Tl)*Cl_MEA_true*1e-6
    Ha = (k2*Cl_MEA_true*Dl_CO2)**.5/kl_CO2

    # KH = E * kl_CO2 / kv_CO2 / (E * kl_CO2 / kv_CO2 + H_CO2_mix)
    KH = Ha * kl_CO2 / kv_CO2 / (Ha * kl_CO2 / kv_CO2 + H_CO2_mix)

    return kv_CO2, kv_H2O, kv_T, KH, [kl_CO2, kv_CO2, kv_H2O, KH, E, Ha], uv, a_e, [ul, uv, h_L, a_e]
