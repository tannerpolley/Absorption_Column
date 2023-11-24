from numpy import polyval
import numpy as np
import pandas as pd
import scipy.linalg as sp
from scipy.stats import linregress
import matplotlib.pyplot as plt


def f_CO2_guess(y_CO2):
    p = [3.58186260e+07, -7.43273232e+06, 5.52655474e+05, -1.79252948e+04, 3.12753126e+02]
    p_linear = [-890.41996267, 138.72207665]
    p_quad = [5216.5, -1535.8, 157.69]

    if y_CO2 < .06:
        sol = round(polyval(p, y_CO2) / 100, 3)
    elif .06 <= y_CO2 <= .08:
        sol = round(polyval(p_linear, y_CO2) / 100, 3)
    elif y_CO2 > .08:
        sol = round(polyval(p_quad, y_CO2) / 100, 3)
    else:
        print('Out of range or wrong type')
        sol = None

    return sol


def fig(i):
    stuff = np.array([93.99, 88.0, 89.55, 90.43, 94.83, 80.49, 62.07, 83.64, 77.5, 93.51, 93.67, 94.73,
 84.7, 86.82, 95.17, 94.39, 92.37])
    return stuff[i]



def secondordfunc(MEAout_xt):
    MEAout_xt2o=np.zeros((MEAout_xt.shape[0],0),dtype=float)
    for ii2 in np.arange(MEAout_xt.shape[1]-1):
        #idxtemp=np.arange(MEAout_xt.shape[1]-1)+ii2*(MEAout_xt.shape[1]-1)
        MEAtempmat=MEAout_xt[:,1:]*MEAout_xt[:,(ii2+1):(ii2+2)]
        MEAout_xt2o=np.concatenate((MEAout_xt2o,MEAtempmat[:,0:(ii2+1)]),axis=1)
    return MEAout_xt2o

def poly_surr(modout_xt, modout_ysim, beta_thresh=0.04):
    # standardize between 0-1
    xt_min = modout_xt.min(axis=0);
    xt_ptp = modout_xt.ptp(axis=0);
    modout_xt = (modout_xt - xt_min) / xt_ptp;
    modout_xt = np.concatenate((0.5 * np.ones((modout_xt.shape[0], 1), dtype=float), modout_xt), axis=1)
    modout_ysim = modout_ysim / 100;
    # first order surrogate

    a = modout_xt.T @ modout_xt
    np.savetxt('data/crap.txt', a)
    a = np.loadtxt('../data/crap.txt')

    beta = sp.inv(a) @ (modout_xt.T @ modout_ysim);
    modout_pred = modout_xt @ beta;
    # second order surrogate
    idx1 = np.squeeze(np.abs(beta) > beta_thresh);
    modout_xtsm = modout_xt[:, np.squeeze(np.abs(beta) > beta_thresh)]
    modout_xt2o = secondordfunc(modout_xtsm);
    modout_xt2 = np.concatenate((modout_xt, modout_xt2o), axis=1)

    b = modout_xt2.T @ modout_xt2
    np.savetxt('data/crap2.txt', b)
    b = np.loadtxt('../data/crap2.txt')

    beta2 = sp.inv(b) @ (modout_xt2.T @ modout_ysim);
    return modout_xt2, idx1, beta2, xt_min, xt_ptp

def poly_surr_pred(modout_xt, idx1, beta2, xt_min, xt_ptp):
    modout_xt = (modout_xt - xt_min) / xt_ptp;
    modout_xt = np.concatenate((0.5 * np.ones((modout_xt.shape[0], 1), dtype=float), modout_xt), axis=1)
    modout_xt2o = secondordfunc(modout_xt[:, idx1])
    modout_xt2 = np.concatenate((modout_xt, modout_xt2o), axis=1)
    pred2 = 100 * modout_xt2 @ beta2;
    return pred2.squeeze()

des1 = np.array(pd.read_csv('../data/run_results_5.csv', index_col=0))
modout_xt = des1[:, 0:5];
modout_ysim = des1[:, -2];

modout_xt2, idx1, beta2, xt_min, xt_ptp = poly_surr(modout_xt, modout_ysim)

def surrogate(X):

    modout_xt1 = np.array([X])
    pred2 = poly_surr_pred(modout_xt1, idx1, beta2, xt_min, xt_ptp)

    return pred2


ypred = []
for i in range(len(modout_xt)):
    ypred.append(surrogate(modout_xt[i]))

df = pd.read_csv('../data/run_results_5.csv', index_col=0)
df['yfit'] = ypred

df.to_csv('data/run_results_5.1.csv')









