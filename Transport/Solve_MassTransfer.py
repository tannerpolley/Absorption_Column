from Parameters import a_p, Clp, Cvp, ϵ, g, R
import numpy as np

log = np.log
exp = np.exp


def solve_masstransfer(rho_mass_l, rho_mass_v, mul_mix, muv_mix, sigma, Dl_CO2, Dv_CO2, Dv_H2O, Dv_T, A, Tl, Tv, ul, uv,
                       H_CO2_mix, Cl_MEA_true):

    # Hold Up
    h_L = np.exp(log(11.4474) + 0.6471 * (log(3.185966) + log(ul) + 1 / 3 * (log(mul_mix) - log(rho_mass_l))))
    h_V = ϵ - h_L

    if h_V < 0:
        print('Error: Flooding as occurred')
        raise TypeError

    d_h = 4 * ϵ / a_p
    Lp = A * a_p / ϵ

    # Compute interfacial area

    a_e = np.exp(log(a_p) + log(1.42)
                 + .116 * (log(rho_mass_l) - log(sigma)
                           + (1 / 3) * log(g)
                           + (4 / 3) * (log(ul) + log(A) - log(Lp))
                           )
                 )

    def f_kl(Dl):
        kl = exp(log(Clp)
                 + (1 / 6) * log(12)
                 + .5 * (log(ul) + log(Dl) - log(h_L) - log(d_h))
                 )
        return kl

    def f_kv(Dv):
        kv = exp(log(Cvp) - log(R) - log(Tv)
                 - 0.5 * log(h_V)
                 + 0.5 * (log(a_p) - log(d_h))
                 + (2 / 3) * log(Dv)
                 + (1 / 3) * (log(muv_mix) - log(rho_mass_v))
                 + (3 / 4) * (log(uv) + log(rho_mass_v) - log(a_p) - log(muv_mix))
                 )
        return kv

    kl_CO2 = f_kl(Dl_CO2)
    kv_CO2 = f_kv(Dv_CO2)
    kv_H2O = f_kv(Dv_H2O)
    kv_T = f_kv(Dv_T) * (R * Tv)

    k2b = 3.95e10 * exp(-6864/Tl) # Greer -
    k2c = 10**(10.49 - 2200/Tl)/1000 # Hakita - Reaction rate of CO2 in aqueous MEA-AMP solution: Experiment and modeling
    k2d = exp(10 - 904.6/Tl)/1000 # Freguia - Modeling of CO2 Capture by Aqueous Monoethanolamine
    k2e = 3.1732e9 * exp(-4936.6 / Tl) * Cl_MEA_true * 1e-6 # Putta, Svendsen, Knuutila 2017 Eqn. 42

    k2 = k2e
    # print(k2b, k2c, k2d, k2e)
    Ha = (k2 * Cl_MEA_true * Dl_CO2) ** .5 / kl_CO2

    KH = Ha * kl_CO2 / kv_CO2 / (Ha * kl_CO2 / kv_CO2 + H_CO2_mix)

    return kv_CO2, kv_H2O, kv_T, KH, [kl_CO2, kv_CO2, kv_H2O, KH, Ha], uv, a_e, [ul, uv, h_L, a_e]
