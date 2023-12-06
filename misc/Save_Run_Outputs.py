from Parameters import sheetnames, make_dfs_dict, n, z
from numpy import arange, round
import numpy as np
import pandas as pd
from BVP.ABS_Column import abs_column
import xlwings as xw

# Put outputs into dictionary and dataframe


def save_run_outputs(Y, Fl_MEA, Fv_N2, Fv_O2, df_param):

    run_type = 'saving'

    outputs_0 = abs_column(z[0], Y.T[0], Fl_MEA, Fv_N2, Fv_O2, df_param, run_type)

    # Initializes each output array in the shape (n, m) where is the # of relevant properties in a group
    # and puts it into a list of output arrays
    output_dict = {}
    for k in sheetnames:
        output_dict[k] = np.zeros((n, len(outputs_0[k])))

    # Updates each output array and the (i, j) height step (i) for relevant group (j)
    for i in range(n):
        outputs = abs_column(z[i], Y.T[i], Fl_MEA, Fv_N2, Fv_O2, df_param, run_type)

        for k in sheetnames:
            output_dict[k][i] = outputs[k]

    # Converts the Outputs dictionary into a dictionary of dataframes
    dfs_dict = make_dfs_dict(output_dict)

    # Updates each sheet in the Excel file with the new data from the df
    wb = xw.Book('data/Results/Profiles.xlsx', read_only=False)

    for sheetname, df in dfs_dict.items():
        try:
            wb.sheets[sheetname].clear()
        except:
            wb.sheets.add(sheetname)
        wb.sheets[sheetname].range("A1").value = df
    for sheet in wb.sheets:
        if sheet.name not in sheetnames:
            sheet.delete()
    wb.save(path=r'data/Results/Profiles.xlsx')
