#' Fit AIPW ComBat model for batch effect correction
#'
#' This function fits an Augmented Inverse Probability Weighting (AIPW) ComBat model 
#' for batch effect correction using an empirical Bayes framework.
#'
#' Note: This function is experimental, and has not been tested on real data. It has only been tested with simulated data with binary (0 or 1) exposures.
#'
#' @param Ys an \code{[n, d]} matrix, for the outcome variables with \code{n} samples in \code{d} dimensions.
#' @param Ts \code{[n]} the labels of the samples, with at most two unique batches.
#' @param Xs \code{[n, r]} the \code{r} covariates/confounding variables, for each of the \code{n} samples, as a data frame with named columns.
#' @param aipw.form A covariate model, given as a formula. Applies for the estimation of propensities for the AIPW step.
#' @param covar.out.form A covariate model, given as a formula. Applies for the outcome regression step of the \code{ComBat} algorithm. Defaults to \code{NULL}, which re-uses \code{aipw.form} for the covariate/outcome model.
#'
#' @return A list containing the fitted AIPW ComBat model and adjusted outcomes:
#' \itemize{
#'    \item{\code{Ys.corrected}} An [n, d] matrix of adjusted outcomes
#'    \item{\code{Model}} A list containing the fitted model components:
#'    \itemize{
#'       \item{\code{Prop_model}} Fitted propensity model
#'       \item{\code{Putcome_models}} List of fitted outcome models for each feature and treatment level
#'       \item{\code{Levels}} Levels of the treatment factor
#'       \item{\code{Treatment_effects}} Estimated treatment effects
#'       \item{\code{AIPW_est}} AIPW estimator results
#'       \item{\code{covar.out.form}} Formula used for the outcome model
#'    }
#' }
#'
#' @import nnet
#'
#' @noRd
cb.learn.fit_aipw_cComBat <- function(Ys, Ts, Xs, aipw.form, covar.out.form) {
  # Estimate propensities
  prop_model <- multinom(Ts ~ ., data = cbind(Xs, Ts = Ts))
  propensities_raw <- predict(prop_model, newdata = Xs, type = "probs")
  
  # Ensure propensities is always a matrix
  if (is.vector(propensities_raw)) {
    propensities <- cbind(1 - propensities_raw, propensities_raw)
    colnames(propensities) <- levels(Ts)
  } else {
    propensities <- propensities_raw
  }
  
  T_levels <- levels(Ts)
  n_features <- ncol(Ys)
  
  outcome_models <- lapply(1:n_features, function(feature) {
    lapply(T_levels, function(t) {
      lm(as.formula(paste("Ys ~", covar.out.form)), 
         data = cbind(data.frame(Ys = Ys[, feature]), Xs)[Ts == t, ])
    })
  })
  
  # Calculate potential outcomes
  potential_outcomes <- array(0, dim = c(nrow(Xs), length(T_levels), n_features))
  for (feature in 1:n_features) {
    for (k in 1:length(T_levels)) {
      potential_outcomes[, k, feature] <- predict(outcome_models[[feature]][[k]], newdata = Xs)
    }
  }
  
  # Calculate AIPW estimator
  n <- nrow(Xs)
  K <- length(T_levels)
  AIPW_est <- array(0, dim = c(K, n_features))
  
  for (k in 1:K) {
    I_k <- as.integer(Ts == T_levels[k])
    for (feature in 1:n_features) {
      AIPW_est[k, feature] <- mean(I_k / propensities[, k] * (Ys[, feature] - potential_outcomes[, k, feature]) + 
                                     potential_outcomes[, k, feature])
    }
  }
  
  # Calculate treatment effects (comparing to the first level)
  treatment_effects <- AIPW_est[-1, , drop = FALSE] - 
    matrix(AIPW_est[1, ], nrow = K-1, ncol = n_features, byrow = TRUE)
  
  # Adjust Y
  Y_star <- matrix(0, nrow = n, ncol = n_features)
  for (feature in 1:n_features) {
    Y_star[, feature] <- Ys[, feature] - 
      potential_outcomes[cbind(1:n, as.integer(Ts), feature)] + 
      potential_outcomes[, 1, feature]
  }
  
  # Prepare estimates
  estimates <- list(
    ATE = treatment_effects,
    potential_outcomes = apply(potential_outcomes, c(2, 3), mean)
  )
  
  # Prepare model object
  model <- list(
    Prop_model = prop_model,
    Outcome_models = outcome_models,
    Levels = T_levels,
    Treatment_effects = treatment_effects,
    AIPW_est = AIPW_est,
    covar.out.form = covar.out.form
  )
  
  return(list(Ys.corrected = Y_star, Model=model))
}


#' Fit AIPW ComBat model for batch effect correction
#'
#' This function applies an Augmented Inverse Probability Weighting (AIPW) ComBat model 
#' for batch effect correction to new data.
#'
#' Note: This function is experimental, and has not been tested on real data. It has only been tested with simulated data with binary (0 or 1) exposures.
#'
#' @param Ys an \code{[n, d]} matrix, for the outcome variables with \code{n} samples in \code{d} dimensions.
#' @param Ts \code{[n]} the labels of the samples, with \code{K < n} levels, as a factor variable.
#' @param Xs \code{[n, r]} the \code{r} covariates/confounding variables, for each of the \code{n} samples, as a data frame with named columns.
#' @param Model a list containing the following parameters:
#'  \itemize{
#'     \item{\code{Prop_model}} Fitted propensity model
#'     \item{\code{Putcome_models}} List of fitted outcome models for each feature and treatment level
#'     \item{\code{Levels}} Levels of the treatment factor
#'     \item{\code{Treatment_effects}} Estimated treatment effects
#'     \item{\code{AIPW_est}} AIPW estimator results
#'     \item{\code{covar.out.form}} Formula used for the outcome model
#'  }
#' This model is output after fitting with \code{\link{cb.correct.aipw_cComBat}}.
#' 
#' @return an \code{[n, d]} matrix, the batch-effect corrected data.
#'
#' @importFrom graphics lines par
#' @importFrom stats cor density dnorm model.matrix pf ppoints prcomp predict
#' qgamma qnorm qqline qqnorm qqplot smooth.spline var
#' @importFrom utils read.delim
#' 
#' @examples
#' library(causalBatch)
#' sim <- cb.sims.sim_linear(a=-1, n=200, err=1/8, unbalancedness=3)
#' # fit batch effect correction for first 100 samples
#' cb.fit <- cb.correct.matching_cComBat(sim$Ys[1:100,,drop=FALSE], sim$Ts[1:100], 
#'                                   data.frame(Covar=sim$Xs[1:100,,drop=FALSE]), "Covar")
#' # apply to all samples
#' cor.dat <- cb.correct.apply_cComBat(sim$Ys, sim$Ts, data.frame(Covar=sim$Xs), cb.fit$Model)
#'
#' @export
cb.correct.apply_aipw_cComBat <- function(Ys, Ts, Xs, Model) {
  # Ensure Ts is a factor
  Ts <- as.factor(Ts)
  # Verify that Ts and the model have the same levels
  for (batch in Ts) {
    if (!(batch %in% Model$Levels)) {
      stop("You have out-of-sample data from batches not in the fit model.")
    }
  }
  
  # Extract model components
  prop_model <- Model$Prop_model
  outcome_models <- Model$Outcome_models
  T_levels <- Model$Levels
  covar.out.form <- Model$covar.out.form
  
  # Estimate propensities for new data
  propensities_raw <- predict(prop_model, newdata = Xs, type = "probs")
  
  # Ensure propensities is always a matrix
  if (is.vector(propensities_raw)) {
    propensities <- cbind(1 - propensities_raw, propensities_raw)
    colnames(propensities) <- T_levels
  } else {
    propensities <- propensities_raw
  }
  
  n_features <- ncol(Ys)
  n <- nrow(Xs)
  
  # Calculate potential outcomes for new data
  potential_outcomes <- array(0, dim = c(nrow(Xs), length(T_levels), n_features))
  for (feature in 1:n_features) {
    for (k in 1:length(T_levels)) {
      potential_outcomes[, k, feature] <- predict(outcome_models[[feature]][[k]], newdata = Xs)
    }
  }
  
  # Adjust Y for new data
  Y_star <- matrix(0, nrow = n, ncol = n_features)
  for (feature in 1:n_features) {
    Y_star[, feature] <- Ys[, feature] - 
      potential_outcomes[cbind(1:n, as.integer(Ts), feature)] + 
      potential_outcomes[, 1, feature]
  }
  
  return(Y_star)
}
