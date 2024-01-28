require(causalBatch)
require(reticulate)

reticulate::use_virtualenv("~/.virtualenvs/neuroharm/", required=TRUE)
reticulate::py_config()
neuroharm <- reticulate::import("neuroHarmonize")

ComBat.GAM <- function(X, batches, covariates, nh.args=NULL) {
  covars <- cbind(data.frame(SITE=batches), covariates)
  res_cgam <- do.call(neuroharm$harmonizationLearn, c(list(X, covars), nh.args))
  return(list(Data=res_cgam[[2]],
              Batches=batches,
              Covariates=covariates,
              Model=res_cgam[[1]]))
}

causal.ComBat.GAM <- function(X, batches, covariates, nh.args=NULL, match.args=list(method="nearest", exact=NULL, replace=FALSE, caliper=.1)) {
  retain.ids <- unique(balance.batches(batches, covariates, match.form, exact=exact, smooth="Age"))
  X.tilde <- X[retain.ids,]; Y.tilde <- covariates[retain.ids,]; t.tilde <- batches[retain.ids]
  
  covars <- cbind(data.frame(SITE=t.tilde), Y.tilde)
  res_cgam <- do.call(neuroharm$harmonizationLearn, c(list(X.tilde, covars), nh.args))
  return(list(Data=res_cgam[[2]],
              Batches=t.tilde,
              Covariates=Y.tilde,
              Retained.Ids=retain.ids,
              Model=res_cgam[[1]]))
}

ComBat.GAM.apply <- function(X, batches, covariates, model) {
  covars <- cbind(data.frame(SITE=batches), covariates)
  
  res_cgam <- do.call(neuroharm$harmonizationApply, list(X, covars, model))
  
  return(list(Data=res_cgam,
              Batches=batches,
              Covariates=covariates,
              Model=model))
}