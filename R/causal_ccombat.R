#' Causal Conditional ComBat
#' 
#' A function for implementing the causal conditional ComBat (causal cComBat) algorithm.
#' This algorithm allows users to remove batch effects (in each dimension), while adjusting for known confounding
#' variables. It is imperative that this function is used in conjunction with domain
#' expertise (e.g., to ensure that the covariates are not colliders, and that the system satisfies the strong
#' ignorability condiiton) to derive causal conclusions. See citation for more details as to the conditions
#' under which conclusions derived are causal.
#' 
#' @importFrom sva ComBat
#' @importFrom stats model.matrix
#' @importFrom stats as.formula
#' @param Ys an \code{[n, d]} matrix, for the outcome variables with \code{n} samples in \code{d} dimensions.
#' @param Ts \code{[n]} the labels of the samples, with \code{K < n} levels, as a factor variable.
#' @param Xs \code{[n, r]} the \code{r} covariates/confounding variables, for each of the \code{n} samples, as a data frame with named columns.
#' @param match.form A formula of columns from \code{Xs}, to be passed directly to \code{\link[MatchIt]{matchit}} for subsequent matching. See \code{formula} argument from \code{\link[MatchIt]{matchit}} for details.
#' @param reference the name of the reference/control batch, against which to match. Defaults to \code{NULL}, which treats the reference batch as the smallest batch.
#' @param match.args A named list arguments for the \code{\link[MatchIt]{matchit}} function, to be used to specify specific matching strategies, where the list names are arguments and the corresponding values the value to be passed to \code{matchit}. Defaults to inexact nearest-neighbor caliper (width 0.1) matching without replacement.
#' @param retain.ratio If the number of samples retained is less than \code{retain.ratio*n}, throws a warning. Defaults to \code{0.05}.
#' @param apply.oos whether to apply batch effect correction to all samples within range of covariate overlap across the different batches. Defaults to \code{FALSE}. 
#' @return a list, containing the following:
#' \itemize{
#'    \item{\code{Ys.corrected}} an \code{[m, d]} matrix, for the \code{m} retained samples in \code{d} dimensions, after correction.
#'    \item{\code{Ts}} \code{[m]} the labels of the \code{m} retained samples, with \code{K < n} levels.
#'    \item{\code{Xs}} the \code{r} covariates/confounding variables for each of the \code{m} retained samples.
#'    \item{\code{Model}} the fit batch effect correction model. See \code{\link[sva]{ComBat}} for details.
#'    \item{\code{InSample.Ids}} the ids which were used to fit the batch effect correction model.
#'    \item{\code{Corrected.Ids}} the ids to which batch effect correction was applied. Differs from \code{InSample.Ids} if \code{apply.oos} is \code{TRUE}.
#' }
#' @param apply.oos A boolean that indicates whether or not to apply the learned batch effect correction to non-matched samples that are still within a region of covariate support. Defaults to \code{FALSE}.
#' 
#' @author Eric W. Bridgeford
#' @references Eric W. Bridgeford, et al. "A Causal Perspective for Batch Effects: When is no answer better than a wrong answer?" Biorxiv (2024). 
#' @references Daniel E. Ho, et al. "MatchIt: Nonparametric Preprocessing for Parametric Causal Inference" JSS (2011). 
#' @references W Evan Johnson, et al. "Adjusting batch effects in microarray expression data using empirical Bayes methods" Biostatistics (2007). 
#' 
#' @section Details:
#' For more details see the help vignette:
#' \code{vignette("causal_ccombat", package = "causalBatch")}
#' 
#' @examples
#' library(causalBatch)
#' sim <- cb.sims.sim_linear(a=-1, n=100, err=1/8, unbalancedness=3)
#' cb.correct.caus_cComBat(sim$Ys, sim$Ts, data.frame(Covar=sim$Xs), "Covar")
#' 
#' @export
cb.correct.caus_cComBat <- function(Ys, Ts, Xs, match.form, reference=NULL, match.args=list(method="nearest", exact=NULL, replace=FALSE, caliper=.1),
                                    retain.ratio=0.05, apply.oos=FALSE) {
  match_obj <- do.call(cb.align.kway_match, list(Ts, Xs, match.form, reference=reference, 
                                                         match.args=match.args, retain.ratio=retain.ratio))
  is.ids <- unique(match_obj$Retained.Ids)
  
  Y.tilde <- Ys[is.ids,,drop=FALSE]; X.tilde <- Xs[is.ids,,drop=FALSE]; T.tilde <- Ts[is.ids]
  
  mod <- model.matrix(as.formula(sprintf("~%s", match.form)), data=X.tilde)
  fit_obj <- causalBatch:::cb.learn.fit_cComBat(Y.tilde, T.tilde, mod = mod)
  fit_obj$Model$Covar.Mod <- match.form
  fit_obj$Model$Reference <- match_obj$Reference
  
  if (apply.oos) {
    oos.ids <- unique(do.call(cb.align.vm_trim, list(Ts, Xs, retain.ratio=retain.ratio)))
    retain.ids <- unique(c(is.ids, oos.ids))
    dat.norm <- cb.correct.apply_cComBat(Ys[retain.ids,,drop=FALSE], Ts[retain.ids], Xs[retain.ids,,drop=FALSE], 
                                         fit_obj$Model)
  } else {
    dat.norm <- fit_obj$Corrected
    retain.ids <- is.ids
  }
  return(list(Ys.corrected=dat.norm,
              Ts=Ts[retain.ids],
              Xs=Xs[retain.ids,,drop=FALSE],
              Model=fit_obj$Model,
              Reference=match_obj$Reference,
              InSample.Ids=is.ids,
              Corrected.Ids=retain.ids))
}

#' K-Way matching
#' 
#' A function for performing k-way matching using the matchIt package. Looks for samples which have corresponding matches across all other treatment levels.
#' 
#' @importFrom magrittr %>%
#' @param Ts \code{[n]} the labels of the samples, with \code{K < n} levels, as a factor variable.
#' @param Xs \code{[n, r]} the \code{r} covariates/confounding variables, for each of the \code{n} samples, as a data frame with named columns.
#' @param match.form A formula of columns from \code{Xs}, to be passed directly to \code{\link[MatchIt]{matchit}} for subsequent matching. See \code{formula} argument from \code{\link[MatchIt]{matchit}} for details.
#' @param reference the name of the reference/control batch, against which to match. Defaults to \code{NULL}, which treats the reference batch as the smallest batch.
#' @param match.args A named list arguments for the \code{\link[MatchIt]{matchit}} function, to be used to specify specific matching strategies, where the list names are arguments and the corresponding values the value to be passed to \code{matchit}. Defaults to inexact nearest-neighbor caliper (width 0.1) matching without replacement.
#' @param retain.ratio If the number of samples retained is less than \code{retain.ratio*n}, throws a warning. Defaults to \code{0.05}.
#' @return a list, containing the following:
#' \itemize{
#'    \item{\code{Retained.Ids}} \code{[m]} vector consisting of the sample ids of the \code{n} original samples that were retained after matching.
#'    \item{\code{Reference}} the reference batch.
#' }
#' 
#' @section Details:
#' For more details see the help vignette:
#' \code{vignette("causal_balancing", package = "causalBatch")}
#' @author Eric W. Bridgeford
#' 
#' @references Eric W. Bridgeford, et al. "A Causal Perspective for Batch Effects: When is no answer better than a wrong answer?" Biorxiv (2024). 
#' @references Daniel E. Ho, et al. "MatchIt: Nonparametric Preprocessing for Parametric Causal Inference" JSS (2011). 
#' 
#' @examples
#' library(causalBatch)
#' sim <- cb.sims.sim_linear(a=-1, n=100, err=1/8, unbalancedness=1.5)
#' cb.align.kway_match(sim$Ts, data.frame(Covar=sim$Xs), "Covar")
#' @export
cb.align.kway_match <- function(Ts, Xs, match.form, reference=NULL, match.args=list(method="nearest", exact=NULL, replace=FALSE, caliper=.1), retain.ratio=0.05) {
  # obtain the smallest batch
  Ts <- as.character(Ts)
  Xs <- cbind(data.frame(Batch=Ts), Xs)
  batch.sum <- Ts %>% table()
  batch.names <- names(batch.sum)
  if (is.null(reference)) {
    tx.batch <- batch.names[which.min(batch.sum)]
  } else {
    tx.batch <- reference
  }
  rownames(Xs) <- 1:nrow(Xs)
  covar.tx <- Xs[Ts == tx.batch,,drop=FALSE]
  
  paired.matches <- lapply(batch.names[batch.names != tx.batch], function(batch) {
    covar.cont <- Xs[Ts == batch,]
    covariate.match(covar.tx, covar.cont, match.form=match.form, match.args=match.args)
  })
  
  # retain control samples with a match
  I.mat <- which(apply(sapply(paired.matches, function(x) x$I.mat.k), c(1), sum) > 0)
  # retain treatment samples with a match
  M.mat <- apply(
    do.call(cbind, lapply(paired.matches, function(x) x$M.mat.k))[I.mat,],
    c(2), sum)
  retain.ids <- as.numeric(c(names(I.mat), names(which(M.mat != 0))))
  
  if (length(retain.ids) < retain.ratio*length(Ts)) {
    warning("Few samples retained by vector matching.")
  }
  
  if (length(retain.ids) == 0) {
    stop("No samples retained by vector matching.")
  }
  
  return(list(Retained.Ids=retain.ids, Reference=tx.batch))
}

#' Pairwise covariate matching
#' 
#' 
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate
#' @importFrom MatchIt matchit
#' @importFrom stats formula
#' @param covar.tx the treatment covariate/label matrix.
#' @param covar.cont the control covariate/label matrix.
#' @param match.form A formula of columns from \code{Xs}, to be passed directly to \code{\link[MatchIt]{matchit}} for subsequent matching. See \code{formula} argument from \code{\link[MatchIt]{matchit}} for details.
#' @param match.args A named list arguments for the \code{\link[MatchIt]{matchit}} function, to be used to specify specific matching strategies, where the list names are arguments and the corresponding values the value to be passed to \code{matchit}.
#' @return a list containing:
#' \itemize{
#'    \item{\code{I.mat.k}} index matrix of control samples that have a match.
#'    \item{\code{M.mat.k}} match matrix of treatment samples that are correspondingly matched to a control.
#' }
#' @noRd
covariate.match <- function(covar.tx, covar.cont, match.form, match.args=NULL) {
  n.kprime <- dim(covar.tx)[1]; n.k <- dim(covar.cont)[1]
  n.matches <- max(1, floor(n.k/n.kprime))
  match <- do.call(matchit, c(list(formula(sprintf("as.factor(Treatment) ~ %s", match.form)),
                                   data=rbind(covar.tx %>% mutate(Treatment = 1),
                                              covar.cont %>% mutate(Treatment = 0)),
                                   ratio=n.matches), match.args))
  mat.mtx <- match$match.matrix
  I.mat.k <- as.numeric(apply(mat.mtx, c(1), function(x) {sum(!is.na(x))}) > 0)
  names(I.mat.k) <- names(match$weights[1:n.kprime])
  M.mat.k <- matrix(FALSE, nrow=n.kprime, ncol=n.k)
  ctrl.names <- names(match$weights[(n.kprime + 1):(n.kprime + n.k)])
  for (i in 1:n.k) {
    for (j in 1:n.kprime) {
      M.mat.k[j,i] <- ifelse(ctrl.names[i] %in% mat.mtx[j,], 1, 0)
    }
  }
  colnames(M.mat.k) <- as.numeric(ctrl.names)
  
  return(list(I.mat.k=I.mat.k, M.mat.k=M.mat.k))
}