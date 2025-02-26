#' Causal Conditional Distance Correlation
#' 
#' A function for implementing the causal conditional distance correlation (causal cDCorr) algorithm.
#' This algorithm allows users to identify whether a treatment causes changes in an outcome, given assorted
#' covariates/confounding variables. It is imperative that this function is used in conjunction with domain
#' expertise (e.g., to ensure that the covariates are not colliders, and that the system satisfies the strong
#' ignorability condiiton) to derive causal conclusions. See citation for more details as to the conditions
#' under which conclusions derived are causal.
#' 
#' @importFrom cdcsis cdcov.test cdcor
#' @importFrom stats as.dist
#' @importFrom stats dist
#' @importFrom stats var
#' @importFrom np npudensbw npseed
#' @param Ys Either:
#' \itemize{
#'    \item{\code{[n, d]} matrix} the outcome variables with \code{n} samples in \code{d} dimensions. In this case, \code{distance} should be \code{FALSE}.
#'    \item{\code{[n, n]} \code{dist} object} a distance object for the \code{n} samples. In this case, \code{distance} should be \code{TRUE}.
#' }
#' @param Ts \code{[n]} the labels of the samples, with \code{K < n} levels, as a factor variable.
#' @param Xs \code{[n, r]} the \code{r} covariates/confounding variables, for each of the \code{n} samples. Ensure that all covariates are encoded properly as continuous or factors (categorical or ordinal) for the bandwidth selection step for kernel width selection.
#' @param prop.form a formula specifying a propensity scoring model. Defaults o \code{NULL}, which uses all of the covariates as exposure predictors. 
#' @param R the number of repetitions for permutation testing. Defaults to \code{1000}.
#' @param dist.method the method used for computing distance matrices. Defaults to \code{"euclidean"}. Other options
#' can be identified by seeing the appropriate documention for the \code{method} argument for the \code{\link[stats]{dist}} function.
#' @param distance a boolean for whether (or not) \code{Ys} are already distance matrices. Defaults to \code{FALSE}, which
#' will use \code{dist.method} parameter to compute an \code{[n, n]} pairwise distance matrix for \code{Ys}.
#' @param width Either:
#' \itemize{
#'  \item{"scott"}{Specifies to use Scott's rule for bandwidth computation.}
#'  \item{"xv"}{Computes bandwidths using cross validation via the \code{\link[np]{npudensbw}} function.}
#'  \item{\code{r} vector}{A vector specifying the bandwidths for each dimension of \code{Xs}.}
#' }
#' Defaults to \code{"scott"}.
#' @param seed a random seed to set. Defaults to \code{1}.
#' @param num.threads The number of threads for parallel processing (if desired). Defaults to \code{1}.
#' @param retain.ratio If the number of samples retained is less than \code{retain.ratio*n}, throws a warning. Defaults to \code{0.05}.
#' @param ddx whether to show additional diagnosis messages. Defaults to \code{FALSE}. Can help with debugging if unexpected results are obtained.
#' @param normalize whether or not to compute the distance correlation (\code{TRUE}) or distance covariance (\code{FALSE}). Skips unnecessary computation if the test outcome is the only item of interest. Defaults to \code{TRUE}.
#' 
#' @return a list, containing the following:
#' \itemize{
#'    \item{\code{Test}} The outcome of the statistical test, from \code{\link[cdcsis]{cdcov.test}}.
#'    \item{\code{Retained.Ids}} The sample indices retained after vertex matching, which correspond to the samples for which statistical inference is performed.
#' }
#' @references Eric W. Bridgeford, et al. "A Causal Perspective for Batch Effects: When is no answer better than a wrong answer?" Imaging Neuroscience (2025). 
#' @references Eric W. Bridgeford, et al. "Learning sources of variability from high-dimensional observational studies" arXiv (2023). 
#' @references Xueqin Wang, et al. "Conditional Distance Correlation" American Statistical Association (2015).
#' 
#' @section Details:
#' For more details see the help vignette:
#' \code{vignette("causal_cdcorr", package = "causalBatch")}
#' 
#' @author Eric W. Bridgeford
#' 
#' @examples
#' library(causalBatch)
#' sim <- cb.sims.sim_linear(a=-1, n=100, err=1/8, unbalancedness=3)
#' cb.detect.caus_cdcorr(sim$Ys, sim$Ts, sim$Xs)
#' 
#' @export
cb.detect.caus_cdcorr <- function(Ys, Ts, Xs, prop.form=NULL, R=1000, dist.method="euclidean",
                                  distance = FALSE, width="scott", seed=1, num.threads=1, retain.ratio=0.05, ddx=FALSE, 
                                  normalize=TRUE) {
  Xs <- as.data.frame(Xs)
  
  # vector match for propensity trimming, and then reduce sub-sample to the
  # propensity matched subset
  retain.ids <- cb.align.vm_trim(Ts, Xs, prop.form=prop.form, retain.ratio=retain.ratio, ddx=ddx)
  if (length(retain.ids) == 0) {
    stop("No samples remain after balancing.")
  }
  if (isTRUE(distance)) {
    DY.tilde <- as.dist(as.matrix(Ys)[retain.ids, retain.ids])
  } else {
    Y.tilde <- Ys[retain.ids,,drop=FALSE]
    DY.tilde = dist(Y.tilde, method=dist.method)
  }
  X.tilde <- Xs[retain.ids,,drop=FALSE]
  if (width == "scott") {
      width <- apply(X.tilde, 2, scotts_rule)
  } else if (width == "xv") {
    if (!is.null(seed)) {
      npseed(seed)
    }
    width <- npudensbw(dat=X.tilde, bwmethod="cv.ml")$bw
  } else if (length(width) != dim(Xs)[2]) {
    stop("You have either not specified a valid bandwidth option. Options are either 'scott', 'xv', or a length-r vector, where r is the number of columns in `Xs`.")
  }
  
  # remove covariate columns with no variance after discarding imbalanced samples
  dropped.cols <- apply(X.tilde, 2, var) == 0
  if (sum(dropped.cols) > 0) {
    message(sprintf("dropping %d columns...", sum(dropped.cols)))
  }
  X.tilde <- data.matrix(X.tilde[, !dropped.cols, drop=FALSE])
  width <- width[!dropped.cols]
  
  DT.tilde <- zero_one_dist(Ts[retain.ids])$DT
  
  # run statistical test
  test.out <- cdcov.test(DY.tilde, DT.tilde, X.tilde, num.bootstrap = R, width = width,
                         seed=seed, num.threads=num.threads, distance=TRUE)
  if (normalize) {
    test.out$statistic <- cdcor(DY.tilde, DT.tilde, X.tilde, width=width, distance=TRUE)$statistic
  }
  return(list(Test=test.out,
              Retained.Ids=retain.ids))
}

#' Vector Matching
#' 
#' A function for implementing the vector matching procedure, a pre-processing step for
#' causal conditional distance correlation. Uses propensity scores to strategically include/exclude
#' samples from subsequent inference, based on whether (or not) there are samples with similar propensity scores
#' across all treatment levels (conceptually, a k-way "propensity trimming"). 
#' It is imperative that this function is used in conjunction with domain expertise to ensure that the covariates are not colliders, 
#' and that the system satisfies the strong ignorability condiiton to derive causal conclusions.
#' 
#' @importFrom nnet multinom
#' @importFrom stats predict
#' 
#' @param Ts \code{[n]} the labels of the samples, with \code{K < n} levels, as a factor variable.
#' @param Xs \code{[n, r]} the \code{r} covariates/confounding variables, for each of the \code{n} samples.
#' @param prop.form a formula specifying a propensity scoring model. Defaults o \code{NULL}, which uses all of the covariates as exposure predictors. 
#' @param retain.ratio If the number of samples retained is less than \code{retain.ratio*n}, throws a warning. Defaults to \code{0.05}.
#' @param ddx whether to show additional diagnosis messages. Defaults to \code{FALSE}. Can help with debugging if unexpected results are obtained.
#' @param reference the name of a reference label, against which to align other labels. Defaults to \code{NULL}, which identifies a shared region of covariate overlap across all labels..
#' @return a \code{[m]} vector containing the indices of samples retained after vector matching.
#' 
#' @references Michael J. Lopez, et al. "Estimation of Causal Effects with Multiple Treatments" Statistical Science (2017). 
#' ran
#' @author Eric W. Bridgeford
#' @section Details:
#' For more details see the help vignette:
#' \code{vignette("causal_balancing", package = "causalBatch")}
#' 
#' @examples
#' library(causalBatch)
#' sim <- cb.sims.sim_linear(a=-1, n=100, err=1/8, unbalancedness=3)
#' cb.align.vm_trim(sim$Ts, sim$Xs)
#' 
#' @export
cb.align.vm_trim <- function(Ts, Xs, prop.form=NULL, retain.ratio=0.05, ddx=FALSE, reference=NULL) {
  Xs = as.data.frame(Xs)
  
  # Fitting the Multinomial Logistic Regression Model
  K <- length(unique(Ts))
  Ts <- as.numeric(factor(Ts, levels=unique(Ts)))
  Ts_unique = unique(Ts)
  
  # Create the right-hand side of the formula based on prop_formula or default to all variables
  if (is.null(prop.form)) {
    rhs_formula <- paste(names(Xs), collapse = " + ")
  } else {
    rhs_formula <- paste(prop.form)
  }
  
  # Create the full formula
  full_formula <- as.formula(paste("factor(Ts) ~", rhs_formula))
  
  # Combine Ts and Xs for the model fitting
  model_data <- cbind(Ts = Ts, Xs)
  
  m <- multinom(full_formula, data = model_data, trace=ddx)
  
  # Making predictions using the fitted model
  pred <- predict(m, newdata = Xs, type = "probs")
  
  # if only binary treatment levels, add a column for the reference
  if (K == 2) {
    pred <- cbind(1 - pred, pred)
    colnames(pred) <- Ts_unique
  }
  
  # Function to calculate the range of predicted probabilities for each treatment
  calculate_Rtable <- function(Tval) {
    Rtab = matrix(0, nrow=2, ncol=length(Ts_unique))
    for (Tp in Ts_unique) {
      Rtab[, Tp] = c(min(pred[Ts == Tp, Tval]), max(pred[Ts == Tp, Tval]))
    }
    c(max(Rtab[1,]), min(Rtab[2,]))
  }
  
  # Creating the Rtable to store the range of predicted probabilities for each treatment
  Rtable <- t(sapply(Ts_unique, calculate_Rtable))
  rownames(Rtable) = Ts_unique
  
  # Function to check if each observation satisfies balance condition for each treatment
  check_balance <- function(i) {
    sapply(as.character(Ts_unique), function(Tval) pred[i, Tval] >= Rtable[Tval, 1] & pred[i, Tval] <= Rtable[Tval, 2])
  }
  
  # Creating the balance_check matrix to check if each observation satisfies balance condition for each treatment
  balance_check <- t(sapply(1:nrow(Xs), check_balance))
  
  # Finding observations that satisfy balance condition for all treatments
  balanced_ids <- apply(balance_check, 1, all)
  retain.ids <- which(balanced_ids)
  
  if (length(retain.ids) < retain.ratio*length(Ts)) {
    warning("Few samples retained by vector matching.")
  }
  
  if (length(retain.ids) == 0) {
    stop("No samples retained by vector matching.")
  }
  
  return(retain.ids)
}

#' A utility to one-hot encode a treatment vector.
#' 
#' @importFrom stats model.matrix
#' @param Ts \code{[n]} the labels of the samples, with \code{K < n} levels, as a factor variable.
#' @param levels optional argument for whether to use pre-specified levels. Defaults to \code{NULL}.
#' @return a list containing the following:
#' \itemize{
#'    \item{\code{ohe}} \code{[n, K]} a one-hot encoding of \code{Ts}.
#'    \item{\code{Levels}} the levels for the one-hot encoding.
#' }
#' @author Eric W. Bridgeford
#' @noRd
ohe <- function(Ts, levels=NULL) {
  # Convert input vector to a factor
  if (!is.null(levels)) {
    Ts_fact <- factor(Ts, levels=levels)
  } else {
    Ts_fact <- factor(Ts)
  }
  
  # Perform one-hot encoding using model.matrix
  encoded_matrix <- model.matrix(~Ts_fact - 1)
  
  # Convert the matrix to a data frame (optional)
  Ts_ohe <- as.data.frame(encoded_matrix)
  return(list(ohe=Ts_ohe, Levels=levels(Ts_fact)))
}

#' A utility to compute the zero-one distances for a treatment vector.
#' 
#' @importFrom stats dist
#' @param Ts \code{[n]} the labels of the samples, with \code{K < n} levels, as a factor variable.
#' @param levels optional argument for whether to use pre-specified levels. Defaults to \code{NULL}.
#' @return a list containing the following:
#' \itemize{
#'    \item{\code{DT}} \code{[n, n]} the pairwise zero-one distance matrix.
#'    \item{\code{Levels}} the levels for the one-hot encoding.
#' }
#' @author Eric W. Bridgeford
#' @noRd
zero_one_dist <- function(Ts, levels=NULL) {
  Ts_ohe <- ohe(Ts, levels=levels)
  
  return(list(DT=dist(Ts_ohe$ohe)/sqrt(2), Levels=Ts_ohe$levels))
}

#'' A utility to compute Scott's rule for bandwidth selection.
#'
#' @param x a \code{r} vector.
#' @return the bandwidth selected.
#' @author Eric W. Bridgeford
#' @noRd
scotts_rule <- function(x) {
  1.06 * sd(x) * length(x)^(-1/5)
}