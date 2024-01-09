from BVP.Run_Model import run_model
from Convert_Data.Get_NCCC_Data import get_NCCC_data
import pandas as pd
import numpy as np
import warnings
warnings.filterwarnings("ignore")


CO2_cap_array = []
results_array = []
inputs_array = []

df = pd.read_csv('data/runs_file_SRP_cases.csv', index_col=0)
for i in range(len(df)):
    try:
        CO2_cap, shooter_message = run_model(df, run=i, show_residuals=False, save_run_results=False)
    except:
        pass

    # inputs_array.append(X)
    # CO2_cap_array.append(CO2_cap)
    # results_array.append(shooter_message)

# data = np.column_stack([Tl_0_guess_array, message_array])
# columns = ['Tl_0_guess', 'message']
# df = pd.DataFrame(data, columns=columns)
# df.index.name = 'Runs'
# df.index += 1
# df.to_csv(r'data/Tl_0_guess_analysis_3.csv')



# data = np.column_stack([inputs_array, CO2_cap_array, results_array])
# columns = ['L', 'G', 'alpha', 'w_MEA', 'y_CO2', 'Tl', 'Tv', 'P', 'Beds', 'CO2 CAP%', 'Results']
# df = pd.DataFrame(data, columns=columns)
# df.index.name = 'Runs'
# df.index += 1
# df.to_csv(r'data/runs_result.csv')
