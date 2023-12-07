from BVP.Run_Model import run_model
from Convert_Data.Get_NCCC_Data import get_NCCC_data
import pandas as pd
import numpy as np


CO2_cap_array = []
results_array = []
inputs_array = []

case_num = 18
for i in range(case_num - 1, case_num):

    X = get_NCCC_data(index=i)
    CO2_cap, shooter_message = run_model(X,
                                         run=i,
                                         show_residuals=False,
                                         save_run_results=True)
    inputs_array.append(X)
    CO2_cap_array.append(CO2_cap)
    results_array.append(shooter_message)

data = np.column_stack([inputs_array, CO2_cap_array, results_array])

columns = ['L', 'G', 'alpha', 'w_MEA', 'y_CO2', 'Tl', 'Tv', 'P', 'Beds', 'CO2 CAP%', 'Results']
df = pd.DataFrame(data, columns=columns)
df.index.name = 'Runs'
df.index += 1
df.to_csv(r'data/runs_result.csv')
