import numpy as np
import pandas as pd
import time
import math
from Parameters import H, n
from BVP.Simulate_Abs_Column import simulate_abs_column
from Convert_Data.Convert_NCCC_Data import convert_NCCC_data
from misc.Save_Run_Outputs import save_run_outputs

np.set_printoptions(suppress=True)


def run_model(X, param_dict=None, run=0, show_info=True, show_residuals=False, save_run_results=True):

    if param_dict is None:
        df_param = pd.read_csv(r'data\Property_Parameters_OG.csv')
        df_param.to_csv(r'data\Property_Parameters.csv')

    else:
        df_param = pd.read_csv(r'data/Property_Parameters_OG.csv')
        for k, v in param_dict.items():
            for i in range(len(df_param[k].to_numpy())):
                if not math.isnan(df_param[k][i]):
                    df_param[k][i] = v[i]
        df_param.to_csv(r'data\Property_Parameters.csv')

    # Converts the NCCC data to dataframes containing all inlet information, especially inlet concentrations
    vapor, liquid = convert_NCCC_data(X, df_param)

    # Collect input data from dataframes
    Fl_z = liquid['n']['CO2'], liquid['n']['MEA'], liquid['n']['H2O']
    Fv_0 = vapor['n']['CO2'], vapor['n']['H2O'], vapor['n']['N2'], vapor['n']['O2']

    Tl_z = X[5]
    Tv_0 = X[6]
    P = X[7]
    n_beds = X[8]
    z = np.linspace(0, H*n_beds, n)

    stages = np.linspace(0, H, n)
    # print(stages)

    inputs = [Fl_z, Fv_0, Tl_z, Tv_0, z, P]

    print(f'Run #{run + 1:03d} --- ', end='')

    # Starts the time tracker for the total computation time for one simulation run
    start = time.time()

    # Simulate the Absorption Column from start to finish given the inlet concentrations of the top liquid and bottom vapor streams
    # Also has a try and except statement to catch runs that ended with too many NaN's

    Y, shooter_message = simulate_abs_column(inputs, df_param, show_residuals)

    # Ends the time tracker for the total computation time for one simulation run
    end = time.time()
    total_time = end - start

    Reset_Counter = int(np.loadtxt(r'counters/Run_Counter.txt'))
    Reset_Counter += 1
    np.savetxt(r'counters/Reset_Counter.txt', np.array([Reset_Counter]))

    if shooter_message == 'This will break the model':
        print(f'Time: {total_time} sec - Run Failed due to convergence issues')
        return 0, 'Run Failed due to convergence issues'

    # Collects data from the final integration output

    Fv_CO2_z, Fv_H2O_z = Y[2, -1], Y[3, -1]
    CO2_cap = abs(Fv_0[0] - Fv_CO2_z) / Fv_0[0] * 100

    Fl_CO2_z, Fl_H2O_z, Tl_z_sim = Y[0, -1], Y[1, -1], Y[4, -1]

    # Computes the relative error between the solution that the shooter found to the actual inlet concentration for the relevant liquid species
    CO2_rel_err = abs(Fl_z[0] - Fl_CO2_z) / Fl_z[0] * 100
    H2O_rel_err = abs(Fl_z[2] - Fl_H2O_z) / Fl_z[2] * 100
    Tl_rel_err = abs(Tl_z - Tl_z_sim) / Tl_z * 100

    # Prints out relevant info such as simulation time, relative errors, CO2% captured, if max iterations were reached, and number of Nan's counted
    if show_info:
        print(f'{CO2_cap:.2f}% - Time: {total_time:0>{4}.1f} sec - % Error: CO2 = {CO2_rel_err:0>{5}.2f}%, H2O = {H2O_rel_err:0>{5}.2f}%, Tl = {Tl_rel_err:0>{5}.2f}% - {shooter_message}')

    # Stores output data into text files (concentrations, mole fractions, and temperatures) (can also plot)
    if save_run_results:
        save_run_outputs(Y, Fl_z[1], Fv_0[2], Fv_0[3], P, df_param, z, stages)

    return np.round(CO2_cap, 2), shooter_message
