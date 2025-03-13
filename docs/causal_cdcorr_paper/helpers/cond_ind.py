import pandas as pd
import numpy as np
import statsmodels.api as sm
import statsmodels
import dodiscover as dod
import hyppo
import scipy as sp
import sklearn as sk
from rpy2.robjects.packages import STAP
from rpy2.robjects import numpy2ri
import os
numpy2ri.activate()

# Get the directory where cond_ind.py is located
CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
R_FILE_PATH = os.path.join(CURRENT_DIR, 'test_statistics.R')

# Load the R file once
with open(R_FILE_PATH, 'r') as f:
    R_STRING = f.read()
COND_IND = STAP(R_STRING, "cond_ind")

def cond_manova(Ys, Ts, Xs, **kwargs):
    stat, pval = COND_IND.test_cmanova(Ys, Ts, Xs)
    return float(pval), float(stat)

def cond_permanova(Ys, Ts, Xs, nrep=1000, **kwargs):
    stat, pval = COND_IND.test_permanova(Ys, Ts, Xs, R=nrep)
    return float(pval), float(stat)

def kcit(Ys, T, Xs, nrep=1000, **kwargs):
    stat, pval = COND_IND.test_kcit(Ys, Ts, Xs, R=nrep)
    return float(pval), float(stat)

def rcit(Ys, Ts, Xs, nrep=1000, **kwargs):
    stat, pval = COND_IND.test_rcit(Ys, Ts, Xs, R=nrep)
    return float(pval), float(stat)

def rcot(Ys, Ts, Xs, nrep=1000, **kwargs):
    stat, pval = COND_IND.test_rcot(Ys, Ts, Xs, R=nrep)
    return float(pval), float(stat)

def wgcm(Ys, Ts, Xs, nrep=1000, **kwargs):
    stat, pval = COND_IND.test_wgcm(Ys, Ts, Xs, R=nrep)
    return float(pval), float(stat)

def gcm(Ys, Ts, Xs, nrep=1000, **kwargs):
    stat, pval = COND_IND.test_gcm(Ys, Ts, Xs, R=nrep)
    return float(pval), float(stat)

def kernelcdtest(Ys, Ts, Xs, R=1000, ncores=1, **kwargs):
    df_dict = {'Group': [int(ti) for ti in Ts]}
    yvars = []
    for i in range(0, Ys.shape[1]):
        yvar = 'Ys{:d}'.format(i)
        df_dict[yvar] = Ys[:,i]
        yvars.append(yvar)
    
    xvars = []
    for i in range(0, Xs.shape[1]):
      xvar = 'X{:d}'.format(i)
      df_dict[xvar] = Xs[:,i]
      xvars.append(xvar)
    
    df = pd.DataFrame(df_dict)
    group_col = 'Group'
    stat, pval = dod.cd.KernelCDTest(null_reps=int(R), n_jobs=int(ncores)).test(df, [group_col], yvars, xvars)
    return pval, stat
    
def cond_dcorr(Ys, Ts, Xs, nrep=1000, ncores=1, **kwargs):
    DT = sk.metrics.pairwise_distances(ohe(Ts), metric="l2")
    DY = sk.metrics.pairwise_distances(Ys, metric="l2")
    if len(Xs.shape) == 1:
        Xs = Xs.reshape(-1, 1)
    cdcorr = hyppo.conditional.ConditionalDcorr(compute_distance=None)
    
    KX = cdcorr._compute_kde(Xs)
    stat, pval = cdcorr.test(DY, DT, KX, reps=nrep, workers=ncores)
    return pval, stat

def dcorr(Ys, Ts, Xs, nrep=1000, ncores=1, **kwargs):
    DT = sk.metrics.pairwise_distances(ohe(Ts), metric="l2")
    DY = sk.metrics.pairwise_distances(Ys, metric="l2")
    stat, pval = hyppo.independence.Dcorr(compute_distance=None, use_cov=False).test(DY, DT, reps=nrep, workers=ncores)
    return pval, stat

def causal_cdcorr(Ys, Ts, Xs, nrep=1000, ncores=1, **kwargs):
    stat, pval = hyppo.causal.CausalCDcorr(use_cov=False).test(Ys, Ts, Xs, reps=nrep, workers=ncores)
    return pval, stat
    
def ohe(Ts):
    K = len(np.unique(Ts))
    ohe_dat = np.zeros((len(Ts), K))
    for t in np.unique(Ts):
        ohe_dat[:,t] = (Ts == t).astype(int)
    return ohe_dat