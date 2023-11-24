import matplotlib.pyplot as plt
from Parameters import z, H, D, dz
import time
import pandas as pd
from numpy import round


def plotting(C_dict, Keys, CO2_PER):

    for i, k in enumerate(Keys):
        if k[1] == 'l':
            k_label = k[:3] + '{' + k[3:] + '}'
            plt.figure(i+1)
            plt.plot(z, C_dict[k], label=f'${k_label}$')
            plt.xlabel('Distance from the top (m)')
            plt.ylabel(f'${k_label}$ mol/m3')

        elif k[1] == 'v':
            k_label = k[:3] + '{' + k[3:] + '}'
            plt.figure(i+1)
            plt.plot(z, C_dict[k], label=f'${k_label}$')
            plt.xlabel('Distance from the top (m)')
            plt.ylabel(f'${k_label}$ mol/m3')
            # plt.gca().invert_xaxis()

        plt.title(f'%CO2 = {CO2_PER:.2f}, H = {H} m, D = {D} m, dz = {dz} steps')
        plt.legend()
        plt.show()




