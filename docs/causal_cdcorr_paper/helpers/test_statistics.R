require(cdcsis)
require(causalBatch)
require(GeneralisedCovarianceMeasure)
require(vegan)
library(reticulate)
require(RCIT)
use_virtualenv("causal")

# Define the Python function in R
py_run_string("
import pandas as pd
import dodiscover as dod

def kernelcdtest(Y, T, X, R=1000):
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
    stat, pval = dod.cd.KernelCDTest(null_reps=int(R)).test(df, [group_col], yvars, xvars)
    return pval, stat
")


test.cdcorr <- function(Ys, Ts, Xs, dist.method="euclidean", width=NULL, R=1000, normalize=TRUE, ...) {
  if (is.null(width)) {
    width <- npudensbw(dat=X.tilde, bwmethod="cv.ml")$bw
  }
  
  DY <- dist(data.matrix(Ys), method=dist.method)
  DT <- causalBatch:::zero_one_dist(Ts)$DT
  Xs <- data.matrix(Xs)
  
  test.out <- tryCatch({
    cdcov.test(DY, DT, Xs, num.bootstrap=R, distance=TRUE, width=width)
  }, error=function(e) {
    return(list(Estimate=NaN, p.value=NaN))
  })
  
  if (normalize) {
    test.out$statistic <- cdcor(DY, DY, Xs, width=width, distance=TRUE)$statistic
  }
  return(list(Estimate=test.out$statistic, p.value=test.out$p.value))
}

test.caus_cdcorr <- function(Ys, Ts, Xs, dist.method="euclidean", R=1000, width=NULL, normalize=TRUE, ...) {
  test.out <- tryCatch({
    cb.detect.caus_cdcorr(Ys, Ts, Xs, R=R, dist.method=dist.method, width=width, normalize=normalize)
  }, error=function(e) {
    return(list(Estimate=NaN, p.value=NaN))
  })
  return(list(Estimate=test.out$Test$statistic, p.value=test.out$Test$p.value))
}

test.kcit <- function(Ys, Ts, Xs, R=1000, ...) {
  res <- tryCatch({
    KCIT(data.matrix(Ys), data.matrix(causalBatch:::ohe(Ts)$ohe), data.matrix(Xs), T_BS=R)},
    error=function(e) {
      return(list(Estimate=NaN, p.value=NaN))
    })
  return(list(Estimate=res$stat, p.value=res$pvalue))
}

test.rcit <- function(Ys, Ts, Xs, R=1000, ...) {
  res <- tryCatch({
    RCIT(data.matrix(Ys), data.matrix(causalBatch:::ohe(Ts)$ohe), data.matrix(Xs), approx="perm", nrep=R)},
    error=function(e) {
      return(list(Estimate=NaN, p.value=NaN))
    })
  return(list(Estimate=res$stat, p.value=res$pvalue))
}

test.rcot <- function(Ys, Ts, Xs, R=1000, ...) {
  res <- tryCatch({
    RCoT(data.matrix(Ys), data.matrix(causalBatch:::ohe(Ts)$ohe), data.matrix(Xs), approx="perm", nrep=R)},
    error=function(e) {
      return(list(Estimate=NaN, p.value=NaN))
    })
  return(list(Estimate=res$stat, p.value=res$pvalue))
}

test.gcm <- function(Ys, Ts, Xs, R=1000, ...) {
  res <- tryCatch({
    gcm.test(data.matrix(Ys), data.matrix(causalBatch:::ohe(Ts)$ohe), data.matrix(Xs), regr.meth="xgboost", nsim=R)},
    error=function(e) {
      return(list(Estimate=NaN, p.value=NaN))
    })
  return(list(Estimate=res$test.statistic, p.value=res$p.value))
}

test.manova <- function(Ys, Ts, Xs, ...) {
  Ys <- data.matrix(Ys); Xs <- data.matrix(Xs)
  Tf <- factor(Ts)
  mod_alt <- lm(Ys ~ Xs + Tf + Tf:Xs)
  mod_null <- lm(Ys ~ Xs)
  
  test <- stats:::anova.mlmlist(mod_alt, mod_null)
  return(list(Estimate=test$Pillai[2], p.value=test$`Pr(>F)`[2]))
}

test.permanova <- function(Ys, Ts, Xs, R=1000, ...) {
  Tf <- factor(Ts)
  Xs <- data.matrix(Xs)
  
  mod_alt <- adonis2(Ys ~ Xs + Tf, method="mahalanobis", by="terms", permutations=R)
  
  return(list(Estimate=mod_alt["Tf", "F"], 
              p.value=mod_alt["Tf", "Pr(>F)"]))
}

test.kcd <- function(Ys, Ts, Xs, R = 1000, ...) {
  # Convert R objects to Python-compatible format
  # Y should be a matrix/data.frame, T and X should be vectors
  Ys <- data.matrix(Ys)
  Xs <- data.matrix(Xs)
  # convert Ts so that it is 0s and 1s
  Ts <- as.numeric(factor(Ts, levels=unique(Ts))) - 1
  
  # Call the Python function using reticulate
  result <- py$kernelcdtest(Ys, Ts, Xs, R)
  
  # Return results as a list
  return(list(
    Estimate = result[[2]],
    p.value = result[[1]]
  ))
}

