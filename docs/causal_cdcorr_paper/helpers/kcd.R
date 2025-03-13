library(reticulate)

use_virtualenv("causal")


# Define the Python function in R
py_run_string("
import pandas as pd
import dodiscover as dod

def kernelcdtest(Y, T, X, R=1000, ncores=1):
    df_dict = {'Group': [int(ti) for ti in T]}
    yvars = []
    for i in range(0, Y.shape[1]):
        yvar = 'Y{:d}'.format(i)
        df_dict[yvar] = Y[:,i]
        yvars.append(yvar)
    
    xvars = []
    for i in range(0, X.shape[1]):
      xvar = 'X{:d}'.format(i)
      df_dict[xvar] = X[:,i]
      xvars.append(xvar)
    
    df = pd.DataFrame(df_dict)
    group_col = 'Group'
    stat, pval = dod.cd.KernelCDTest(null_reps=int(R), n_jobs=int(ncores)).test(df, [group_col], yvars, xvars)
    return pval, stat
")