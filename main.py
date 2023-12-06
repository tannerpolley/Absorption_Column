from BVP.Run_Model import run_model
from Convert_Data.Get_NCCC_Data import get_NCCC_data
import pandas as pd
import numpy as np


CO2_cap_array = []
results_array = []
inputs_array = []

# avg = np.array([1.88, .595, .167, .297, .099])
# std = np.array([.9, .2, .1, .02, .015])
# l_bounds = avg - std
# u_bounds = avg + std
# sampler = qmc.LatinHypercube(d=5)
# sample = sampler.random(n=500)
# sample = qmc.scale(sample, l_bounds, u_bounds)
# data = sample

# df = pd.read_csv('data/test_for_failures_outputv2[8].csv')

case_num = 18
for i in range(case_num - 1, case_num):

    # m_T_l, m_T_v, alpha, w_MEA, y_CO2
    X = get_NCCC_data(index=i)
    CO2_cap, shooter_message = run_model(X,
                                         run=i,
                                         show_residuals=False,
                                         save_run_results=True)
    inputs_array.append(X)
    CO2_cap_array.append(CO2_cap)
    results_array.append(shooter_message)

data = np.column_stack([inputs_array, CO2_cap_array, results_array])

columns = ['L', 'G', 'alpha', 'w_MEA', 'y_CO2', 'CO2 CAP%', 'Results']
df = pd.DataFrame(data, columns=columns)
df.index.name = 'Runs'
df.index += 1
df.to_csv(r'data/runs_result.csv')
