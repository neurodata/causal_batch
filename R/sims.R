#' Sigmoidal Simulation
#' 
#' @importFrom stats rbinom
#' @importFrom stats rnorm
#' @param n the number of samples. Defaults to \code{100}.
#' @param pi the balance between the classes, where samples will be from group 1
#' with probability \code{pi}, and group 2 with probability \code{1 - pi}. Defaults
#' to \code{0.5}.
#' @param eff_sz the treatment effect between the different groups. Defaults to \code{1}.
#' @param alpha the alpha for the covariate sampling procedure. Defaults to \code{2}.
#' @param unbalancedness the level of covariate dissimilarity between the covariates
#' for each of the groups. Defaults to \code{1}.
#' @param null whether to generate a null simulation. Defaults to \code{FALSE}. Same behavior can be achieved by setting \code{eff_sz = 0}.
#' @param a the first parameter for the covariate/outcome relationship. Defaults to \code{-4}.
#' @param b the second parameter for the covariate/outcome relationship. Defaults to \code{8}.
#' @param err the level of noise for the simulation. Defaults to \code{1/2}.
#' @param nbreaks the number of breakpoints for computing the expected outcome at a given covariate level
#' for each batch. Defaults to \code{200}.
#' @return a list, containing the following:
#' \item{Y}{an \code{[n, 2]} matrix, containing the outcomes for each sample. The first dimension contains the "treatment effect".}
#' \item{Ts}{an \code{[n, 1]} matrix, containing the group/batch labels for each sample.}
#' \item{Xs}{an \code{[n, 1]} matrix, containing the covariate values for each sample.}
#' \item{Eps}{an \code{[n, 1]} matrix, containing the error for each sample.}
#' \item{x.bounds}{the theoretical bounds for the covariate values.}
#' \item{Ytrue}{an \code{[nbreaks*2, 2]} matrix, containing the expected outcomes at a covariate level indicated by \code{Xtrue}.}
#' \item{Ttrue}{an \code{[nbreaks*2,1]} matrix, indicating the group/batch the expected outcomes and covariate breakpoints correspond to.}
#' \item{Xtrue}{an \code{[nbreaks*2, 1]} matrix, indicating the values of the covariate breakpoints for the theoretical expected outcome in \code{Ytrue}.}
#' \item{Overlap}{the theoretical degree of overlap between the covariate distributions for each of the two groups/batches.}
#' 
#' @section Details:
#' 
#' A sigmoidal relationship between the covariate and the outcome. The first dimension of the outcome is:
#' \deqn{Y_i = a\times \text{sigmoid}(b \times X_i) - a - \text{eff\_sz} \times T_i + \frac{1}{2} \epsilon_i}
#' 
#' where the batch/group labels are:
#' \deqn{T_i \overset{iid}{\sim} Bern(\pi)}
#' 
#' The beta coefficient for the covariate sampling is:
#' \deqn{\beta = \alpha \times \text{unbalancedness}}
#' The covariate values for the first batch are:
#' \deqn{X_i | T_i = 0 \overset{ind}{\sim} 2 Beta(\alpha, \beta) - 1}
#' and the covariate values for the second batch are:
#' \deqn{X_i | T_i = 1 \overset{ind}{\sim} 2 Beta(\beta, \alpha) - 1}
#' Finally, the error terms are:
#' \deqn{\epsilon_i \overset{iid}{\sim} Norm(0, \text{err}^2)}
#' 
#' For more details see the help vignette:
#' \code{vignette("causal_simulations", package = "causalBatch")}
#' 
#' @references Eric W. Bridgeford, et al. "A Causal Perspective for Batch Effects: When is no answer better than a wrong answer?" Biorxiv (2024). 
#' 
#' @author Eric W. Bridgeford
#' 
#' @examples
#' 
#' library(causalBatch)
#' sim = cb.sims.sim_sigmoid()
#' 
#' @export
cb.sims.sim_sigmoid <- function(n=100, pi=.5, eff_sz=1, alpha=2, unbalancedness=1, null=FALSE, 
                        a=-4, b=8, err=1/2, nbreaks=200) {
  batches <- rbinom(n=n, size=1, prob=pi)
  
  beta <- alpha*unbalancedness
  xs <- cb.sims.covar_generator(batches, alpha, beta, beta, alpha)
  eps <- err*rnorm(n=n, mean=0, sd=1)
  ys <- sapply(xs, function(x) {
    a*sigmoid(b*x)
  }) + eps - a
  
  xtrue_tmp <- seq(from=-1, to=1, length.out=nbreaks)
  xtrue <- c(xtrue_tmp, xtrue_tmp)
  batch_true <- c(rep(0, nbreaks), rep(1, nbreaks))
  ytrue <- sapply(xtrue, function(x) {
    a*sigmoid(b*x) - a
  })
  
  if (!null) {
    ys <- ys - eff_sz*batches
    ytrue <- ytrue - eff_sz*batch_true
  }
  
  # ys <- ys * -2*(batches - .5)
  
  return(list(Ys=cbind(matrix(ys, ncol=1), rnorm(n)), Ts=matrix(batches, ncol=1), Xs=matrix(xs, ncol=1),
              Eps=matrix(eps, ncol=1), x.bounds=c(-1, 1), Ytrue=cbind(ytrue, 0), Xtrue=xtrue, 
              Ttrue=batch_true, Overlap=cb.sims.get_beta_overlap(alpha, beta, beta, alpha)))
}

#' Linear Simulation
#' 
#' @importFrom stats rbinom
#' @importFrom stats rnorm
#' @param n the number of samples. Defaults to \code{100}.
#' @param pi the balance between the classes, where samples will be from group 1
#' with probability \code{pi}, and group 2 with probability \code{1 - pi}. Defaults
#' to \code{0.5}.
#' @param eff_sz the treatment effect between the different groups. Defaults to \code{1}.
#' @param alpha the alpha for the covariate sampling procedure. Defaults to \code{2}.
#' @param unbalancedness the level of covariate dissimilarity between the covariates
#' for each of the groups. Defaults to \code{1}.
#' @param null whether to generate a null simulation. Defaults to \code{FALSE}. Same behavior can be achieved by setting \code{eff_sz = 0}.
#' @param a the first parameter for the covariate/outcome relationship. Defaults to \code{-2}.
#' @param b the second parameter for the covariate/outcome relationship. Defaults to \code{-1}.
#' @param err the level of noise for the simulation. Defaults to \code{1/2}.
#' @param nbreaks the number of breakpoints for computing the expected outcome at a given covariate level
#' for each batch. Defaults to \code{200}.
#' @return a list, containing the following:
#' \item{Ys}{an \code{[n, 2]} matrix, containing the outcomes for each sample. The first dimension contains the "treatment effect".}
#' \item{Ts}{an \code{[n, 1]} matrix, containing the group/batch labels for each sample.}
#' \item{Xs}{an \code{[n, 1]} matrix, containing the covariate values for each sample.}
#' \item{Eps}{an \code{[n, 1]} matrix, containing the error for each sample.}
#' \item{x.bounds}{the theoretical bounds for the covariate values.}
#' \item{Ytrue}{an \code{[nbreaks*2, 2]} matrix, containing the expected outcomes at a covariate level indicated by \code{Xtrue}.}
#' \item{Ttrue}{an \code{[nbreaks*2,1]} matrix, indicating the group/batch the expected outcomes and covariate breakpoints correspond to.}
#' \item{Xtrue}{an \code{[nbreaks*2, 1]} matrix, indicating the values of the covariate breakpoints for the theoretical expected outcome in \code{Ytrue}.}
#' \item{Overlap}{the theoretical degree of overlap between the covariate distributions for each of the two groups/batches.}
#' 
#' @section Details:
#' 
#' A linear relationship between the covariate and the outcome. The first dimension of the outcome is:
#' \deqn{Y_i = a\times (X_i + b) - \text{eff\_sz} \times T_i + \frac{1}{2} \epsilon_i}
#' 
#' where the batch/group labels are:
#' \deqn{T_i \overset{iid}{\sim} Bern(\pi)}
#' 
#' The beta coefficient for the covariate sampling is:
#' \deqn{\beta = \alpha \times \text{unbalancedness}}
#' The covariate values for the first batch are:
#' \deqn{X_i | T_i = 0 \overset{ind}{\sim} 2 Beta(\alpha, \beta) - 1}
#' and the covariate values for the second batch are:
#' \deqn{X_i | T_i = 1 \overset{ind}{\sim} 2 Beta(\beta, \alpha) - 1}
#' Finally, the error terms are:
#' \deqn{\epsilon_i \overset{iid}{\sim} Norm(0, \text{err}^2)}
#' 
#' For more details see the help vignette:
#' \code{vignette("causal_simulations", package = "causalBatch")}
#' 
#' @references Eric W. Bridgeford, et al. "A Causal Perspective for Batch Effects: When is no answer better than a wrong answer?" Biorxiv (2024). 
#' 
#' @author Eric W. Bridgeford
#' 
#' @examples
#' 
#' library(causalBatch)
#' sim = cb.sims.sim_linear()
#' 
#' @export
cb.sims.sim_linear <- function(n=100, pi=.5, eff_sz=1, alpha=2, unbalancedness=1, err=1/2, 
                              null=FALSE, a=-2, b=-1, nbreaks=200) {
  batches <- rbinom(n=n, size=1, prob=pi)
  
  beta <- alpha*unbalancedness
  xs <- cb.sims.covar_generator(batches, alpha, beta, beta, alpha)
  
  eps <- err*rnorm(n=n, mean=0, sd=1)
  ys <- a*(xs + b) + eps
  
  xtrue_tmp <- seq(from=-1, to=1, length.out=nbreaks)
  xtrue <- c(xtrue_tmp, xtrue_tmp)
  batch_true <- c(rep(0, nbreaks), rep(1, nbreaks))
  ytrue <- a*(xtrue + b)
  
  if (!null) {
    ys <- ys - eff_sz*batches
    ytrue <- ytrue - eff_sz*batch_true
  }
  
  return(list(Ys=cbind(matrix(ys, ncol=1), rnorm(n)), Ts=matrix(batches, ncol=1), Xs=matrix(xs, ncol=1),
              Eps=matrix(eps, ncol=1), x.bounds=c(-1, 1), Ytrue=cbind(ytrue, 0), Xtrue=xtrue,  
              Ttrue=batch_true, Overlap=cb.sims.get_beta_overlap(alpha, beta, beta, alpha)))
}

#' Impulse Simulation
#' 
#' @importFrom stats rbinom
#' @importFrom stats rnorm
#' @importFrom stats dnorm
#' @param n the number of samples. Defaults to \code{100}.
#' @param pi the balance between the classes, where samples will be from group 1
#' with probability \code{pi}, and group 2 with probability \code{1 - pi}. Defaults
#' to \code{0.5}.
#' @param eff_sz the treatment effect between the different groups. Defaults to \code{1}.
#' @param alpha the alpha for the covariate sampling procedure. Defaults to \code{2}.
#' @param unbalancedness the level of covariate dissimilarity between the covariates
#' for each of the groups. Defaults to \code{1}.
#' @param null whether to generate a null simulation. Defaults to \code{FALSE}. Same behavior can be achieved by setting \code{eff_sz = 0}.
#' @param a the first parameter for the covariate/outcome relationship. Defaults to \code{-0.5}.
#' @param b the second parameter for the covariate/outcome relationship. Defaults to \code{1/2}.
#' @param c the third parameter for the covariate/outcome relationship. Defaults to \code{1}.
#' @param err the level of noise for the simulation. Defaults to \code{1/2}.
#' @param nbreaks the number of breakpoints for computing the expected outcome at a given covariate level
#' for each batch. Defaults to \code{200}.
#' @return a list, containing the following:
#' \item{Ys}{an \code{[n, 2]} matrix, containing the outcomes for each sample. The first dimension contains the "treatment effect".}
#' \item{Ts}{an \code{[n, 1]} matrix, containing the group/batch labels for each sample.}
#' \item{Xs}{an \code{[n, 1]} matrix, containing the covariate values for each sample.}
#' \item{Eps}{an \code{[n, 1]} matrix, containing the error for each sample.}
#' \item{x.bounds}{the theoretical bounds for the covariate values.}
#' \item{Ytrue}{an \code{[nbreaks*2, 2]} matrix, containing the expected outcomes at a covariate level indicated by \code{Xtrue}.}
#' \item{Ttrue}{an \code{[nbreaks*2,1]} matrix, indicating the group/batch the expected outcomes and covariate breakpoints correspond to.}
#' \item{Xtrue}{an \code{[nbreaks*2, 1]} matrix, indicating the values of the covariate breakpoints for the theoretical expected outcome in \code{Ytrue}.}
#' \item{Overlap}{the theoretical degree of overlap between the covariate distributions for each of the two groups/batches.}
#' 
#' @section Details:
#' 
#' A sigmoidal relationship between the covariate and the outcome. The first dimension of the outcome is:
#' \deqn{Y_i = c \times \phi(X_i, \mu = a, \sigma = b) - \text{eff\_sz} \times T_i + \frac{1}{2} \epsilon_i}
#' where \eqn{\phi(x, \mu, \sigma)} is the probability density function for the normal distribution with
#' mean \eqn{\mu} and standard deviation \eqn{\sigma}.
#' 
#' where the batch/group labels are:
#' \deqn{T_i \overset{iid}{\sim} Bern(\pi)}
#' 
#' The beta coefficient for the covariate sampling is:
#' \deqn{\beta = \alpha \times \text{unbalancedness}}
#' The covariate values for the first batch are:
#' \deqn{X_i | T_i = 0 \overset{ind}{\sim} 2 Beta(\alpha, \beta) - 1}
#' and the covariate values for the second batch are:
#' \deqn{X_i | T_i = 1 \overset{ind}{\sim} 2 Beta(\beta, \alpha) - 1}
#' Note that \eqn{X_i | T_i = 0 \overset{D}{=} - X_i | T_i = 1}, or that the covariates are symmetric
#' about the origin in distribution.
#' 
#' Finally, the error terms are:
#' \deqn{\epsilon_i \overset{iid}{\sim} Norm(0, \text{err}^2)}
#' 
#' For more details see the help vignette:
#' \code{vignette("causal_simulations", package = "causalBatch")}
#' 
#' @references Eric W. Bridgeford, et al. "A Causal Perspective for Batch Effects: When is no answer better than a wrong answer?" Biorxiv (2024). 
#' 
#' @author Eric W. Bridgeford
#' 
#' @examples
#' 
#' library(causalBatch)
#' sim = cb.sims.sim_impulse()
#' 
#' @export
cb.sims.sim_impulse <- function(n=100, pi=.5, eff_sz=1, alpha=2, unbalancedness=1,
                                err=1/2, null=FALSE, a=-.5, b=1/2, c=4, nbreaks=200) {
  batches <- rbinom(n=n, size=1, prob=pi)
  
  beta <- alpha*unbalancedness
  xs <- cb.sims.covar_generator(batches, alpha, beta, beta, alpha)
  
  eps <- 1/2*rnorm(n=n, mean=0, sd=1)
  ys <- dnorm(xs, mean=a, sd=b)*c + eps
  
  xtrue_tmp <- seq(from=-1, to=1, length.out=nbreaks)
  xtrue <- c(xtrue_tmp, xtrue_tmp)
  batch_true <- c(rep(0, nbreaks), rep(1, nbreaks))
  ytrue <- dnorm(xtrue, mean=a, sd=b)*c
  
  if (!null) {
    ys <- ys - eff_sz*batches
    ytrue <- ytrue - eff_sz*batch_true
  }
  
  return(list(Ys=cbind(matrix(ys, ncol=1), rnorm(n)), Ts=matrix(batches, ncol=1), Xs=matrix(xs, ncol=1),
              Eps=matrix(eps, ncol=1), x.bounds=c(-1, 1), Ytrue=cbind(ytrue, 0), Xtrue=xtrue,  
              Ttrue=batch_true, Overlap=cb.sims.get_beta_overlap(alpha, beta, beta, alpha)))
}

#' Impulse Simulation with Asymmetric Covariates
#' 
#' @importFrom stats rbinom
#' @importFrom stats rnorm
#' @importFrom stats dnorm
#' @param n the number of samples. Defaults to \code{100}.
#' @param pi the balance between the classes, where samples will be from group 1
#' with probability \code{pi}, and group 2 with probability \code{1 - pi}. Defaults
#' to \code{0.5}.
#' @param eff_sz the treatment effect between the different groups. Defaults to \code{1}.
#' @param alpha the alpha for the covariate sampling procedure. Defaults to \code{2}.
#' @param unbalancedness the level of covariate dissimilarity between the covariates
#' for each of the groups. Defaults to \code{1}.
#' @param null whether to generate a null simulation. Defaults to \code{FALSE}. Same behavior can be achieved by setting \code{eff_sz = 0}.
#' @param a the first parameter for the covariate/outcome relationship. Defaults to \code{-0.5}.
#' @param b the second parameter for the covariate/outcome relationship. Defaults to \code{1/2}.
#' @param c the third parameter for the covariate/outcome relationship. Defaults to \code{1}.
#' @param err the level of noise for the simulation. Defaults to \code{1/2}.
#' @param nbreaks the number of breakpoints for computing the expected outcome at a given covariate level
#' for each batch. Defaults to \code{200}.
#' @return a list, containing the following:
#' \item{Ys}{an \code{[n, 2]} matrix, containing the outcomes for each sample. The first dimension contains the "treatment effect".}
#' \item{Ts}{an \code{[n, 1]} matrix, containing the group/batch labels for each sample.}
#' \item{Xs}{an \code{[n, 1]} matrix, containing the covariate values for each sample.}
#' \item{Eps}{an \code{[n, 1]} matrix, containing the error for each sample.}
#' \item{x.bounds}{the theoretical bounds for the covariate values.}
#' \item{Ytrue}{an \code{[nbreaks*2, 2]} matrix, containing the expected outcomes at a covariate level indicated by \code{Xtrue}.}
#' \item{Ttrue}{an \code{[nbreaks*2,1]} matrix, indicating the group/batch the expected outcomes and covariate breakpoints correspond to.}
#' \item{Xtrue}{an \code{[nbreaks*2, 1]} matrix, indicating the values of the covariate breakpoints for the theoretical expected outcome in \code{Ytrue}.}
#' \item{Overlap}{the theoretical degree of overlap between the covariate distributions for each of the two groups/batches.}
#' 
#' @section Details:
#' 
#' A sigmoidal relationship between the covariate and the outcome. The first dimension of the outcome is:
#' \deqn{Y_i = c \times \phi(X_i, \mu = a, \sigma = b) - \text{eff\_sz} \times T_i + \frac{1}{2} \epsilon_i}
#' where \eqn{\phi(x, \mu, \sigma)} is the probability density function for the normal distribution with
#' mean \eqn{\mu} and standard deviation \eqn{\sigma}.
#' 
#' where the batch/group labels are:
#' \deqn{T_i \overset{iid}{\sim} Bern(\pi)}
#' 
#' The beta coefficient for the covariate sampling is:
#' \deqn{\beta = \alpha \times \text{unbalancedness}}
#' The covariate values for the first batch are asymmetric, in that for the first batch:
#' \deqn{X_i | T_i = 0 \overset{ind}{\sim} 2 Beta(\alpha, \alpha) - 1}
#' and the covariate values for the second batch are:
#' \deqn{X_i | T_i = 1 \overset{ind}{\sim} 2 Beta(\beta, \alpha) - 1}
#' Finally, the error terms are:
#' \deqn{\epsilon_i \overset{iid}{\sim} Norm(0, \text{err}^2)}
#' 
#' For more details see the help vignette:
#' \code{vignette("causal_simulations", package = "causalBatch")}
#' 
#' @references Eric W. Bridgeford, et al. "A Causal Perspective for Batch Effects: When is no answer better than a wrong answer?" Biorxiv (2024). 
#' 
#' @author Eric W. Bridgeford
#' 
#' @examples
#' 
#' library(causalBatch)
#' sim = cb.sims.sim_impulse_asycov()
#' 
#' @export
cb.sims.sim_impulse_asycov <- function(n=100, pi=.5, eff_sz=1, alpha=2, unbalancedness=1, 
                               null=FALSE, a=-.5, b=1/2, c=4, err=1/2, nbreaks=200) {
  batches <- rbinom(n=n, size=1, prob=pi)
  
  beta <- alpha*unbalancedness
  xs <- cb.sims.covar_generator(batches, alpha, alpha, beta, alpha)
  
  eps <- err*rnorm(n=n, mean=0, sd=1)
  ys <- dnorm(xs, mean=a, sd=b)*c + eps
  
  xtrue_tmp <- seq(from=-1, to=1, length.out=nbreaks)
  xtrue <- c(xtrue_tmp, xtrue_tmp)
  batch_true <- c(rep(0, nbreaks), rep(1, nbreaks))
  ytrue <- dnorm(xtrue, mean=a, sd=b)*c
  
  if (!null) {
    ys <- ys - eff_sz*batches
    ytrue <- ytrue - eff_sz*batch_true
  }
  
  return(list(Ys=cbind(matrix(ys, ncol=1), rnorm(n)), Ts=matrix(batches, ncol=1), Xs=matrix(xs, ncol=1),
              Eps=matrix(eps, ncol=1), x.bounds=c(-1, 1), Ytrue=cbind(ytrue, 0), Xtrue=xtrue, 
              Ttrue=batch_true, Overlap=cb.sims.get_beta_overlap(alpha, alpha, beta, alpha)))
}

#' Compute overlap of two beta distributions
#'
#' @importFrom stats dbeta
#' @param a1 alpha of the first covariate distribution.
#' @param b1 beta of the first covariate distribution.
#' @param a2 alpha of the second covariate distribution.
#' @param b2 beta of the second covariate distribution.
#' @param nbreaks the number of breakpoints for approximating the covariate overlap.
#' @return the level of covariate overlap, corresponding to the AUC upper-bounded
#' by the probability density functions for each of the beta distributions.
cb.sims.get_beta_overlap <- function(a1, b1, a2, b2, nbreaks=1000) {
  xbreaks=seq(from=-1, to=1, length.out=nbreaks)
  dbeta1 <- dbeta(1/2*(xbreaks + 1), a1, b1)
  dbeta1 <- dbeta1/sum(dbeta1)
  dbeta2 <- dbeta(1/2*(xbreaks + 1), a2, b2)
  dbeta2 <- dbeta2/sum(dbeta2)
  sum(pmin(dbeta1, dbeta2))
}

#' Sigmoid function
#' 
#' @param x the value
#' @return sigmoid(x)
sigmoid <- function(x) {
  1/(1 + exp(-x))
}

#' Covariate generator function
#'
#' @importFrom stats rbeta 
#' @param batches an \code{n} vector, consisting of the batch labels for each of the \code{n} samples.
#' @param a1 alpha of the first covariate distribution.
#' @param b1 beta of the first covariate distribution.
#' @param a2 alpha of the second covariate distribution.
#' @param b2 beta of the second covariate distribution.
#' @return an \code{n} vector, consisting of the covariate values for each of the \code{n} samples.
cb.sims.covar_generator <- function(batches, a1, b1, a2, b2) {
  2*sapply(batches, function(batch) {
    ifelse(batch==0, rbeta(1, a1, b1), rbeta(1, a2, b2))
  }) - 1
}

