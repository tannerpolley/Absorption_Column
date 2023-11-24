import numpy as np
from Parameters import liquid_species, vapor_species


def heat_capacity(T, phase, x, w):
    if phase == 'vapor':



        coefficients = {'CO2': np.array([19.1, 7.342e-2, -5.602e-5, 1.715e-8]),
                        'H2O': np.array([32.2, 1.924e-3, 1.055e-5, -3.596e-9]),
                        'N2': np.array([31.15, -1.357e-2, 2.68e-5, -1.168e-8]),
                        'O2': np.array([28.11, -3.68e-6, 1.746e-5, -1.065e-8]),
                        }
        Cp = []
        for sp in vapor_species:
            C1, C2, C3, C4 = coefficients[sp]
            Cp.append(C1 + C2 * T + C3 * T ** 2 + C4 * T ** 3)

        return Cp

    elif phase == 'liquid':
        # coefficients = {'CO2': np.array([-203.6060, 1523.290, 3196.413, 2474.455, 3.855326]),
        #                 'MEA': np.array([520.3, 2749, 5313, 0, 0]),
        #                 'H2O': np.array([-203.6060, 1523.290, 3196.413, 2474.455, 3.855326]),
        #                 'N2': np.array([0, 0, 0, 0, 0]),
        #                 'O2': np.array([0, 0, 0, 0, 0]),
        #                 'MEAH': np.array([520.3, 2749, 5313, 0, 0]),
        #                 'MEACOO': np.array([520.3, 2749, 5313, 0, 0]),
        #                 'HCO3': np.array([-203.6060, 1523.290, 3196.413, 2474.455, 3.855326]),
        #                 'OH': np.array([-203.6060, 1523.290, 3196.413, 2474.455, 3.855326]),
        #                 'H3O': np.array([-203.6060, 1523.290, 3196.413, 2474.455, 3.855326])
        #                 }

        coefficients = {'CO2': np.array([276370, -2090.1, 8.125, -.014116, 9.3701e-6]),
                        'MEA': np.array([114000, 158.6, 0, 0, 0]),
                        'H2O': np.array([276370, -2090.1, 8.125, -.014116, 9.3701e-6]),
                        }
        Cp = []
        for sp in liquid_species:
            A, B, C, D, E = coefficients[sp]
            Cp.append((A + B * T + C * T ** 2 + D * T ** 3 + E * T ** 4)/1000)

        Cpl_CO2 = (x[2]*Cp[2] + x[1]*Cp[1])/x[0]*(1/(1 - w[0]) - 1)

        Cp[0] = Cpl_CO2

        return Cp

