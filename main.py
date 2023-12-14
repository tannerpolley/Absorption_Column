from BVP.Run_Model import run_model
from Convert_Data.Get_NCCC_Data import get_NCCC_data
import pandas as pd
import numpy as np


CO2_cap_array = []
results_array = []
inputs_array = []

case_num = 18
Tl_0_guess_array = np.linspace(323.75400801000, 323.75400805000, 5001)
message_array = []
for i, Tl_0_guess in enumerate(Tl_0_guess_array):

    # X = get_NCCC_data(index=i)
    X = [1.89,	0.6308,	0.141,	0.302,	0.1019,	315.39,	319.22,	109180,	1]

    message = run_model(X, Tl_0_guess, run=i, show_residuals=False, save_run_results=False, show_info=False)
    message_array.append(message)

    # CO2_cap, shooter_message = run_model(X,
    #                                      Tl_0_guess,
    #                                      run=i,
    #                                      show_residuals=False,
    #                                      save_run_results=True)
    # inputs_array.append(X)
    # CO2_cap_array.append(CO2_cap)
    # results_array.append(shooter_message)

data = np.column_stack([Tl_0_guess_array, message_array])
columns = ['Tl_0_guess', 'message']
df = pd.DataFrame(data, columns=columns)
df.index.name = 'Runs'
df.index += 1
df.to_csv(r'data/Tl_0_guess_analysis.csv')



# data = np.column_stack([inputs_array, CO2_cap_array, results_array])
# columns = ['L', 'G', 'alpha', 'w_MEA', 'y_CO2', 'Tl', 'Tv', 'P', 'Beds', 'CO2 CAP%', 'Results']
# df = pd.DataFrame(data, columns=columns)
# df.index.name = 'Runs'
# df.index += 1
# df.to_csv(r'data/runs_result.csv')
