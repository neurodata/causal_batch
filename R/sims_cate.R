#' Simulate Binary Covariates
#'
#' @importFrom stats rbinom
#' @param groups group assignments, as an \code{[n]} vector, which take values from \code{0} to \code{K} for \code{K} groups.
#' @param balance a parameter that governs the similarity between the binary covariate distributions for the \code{0} group
#' against the other \code{K-1} groups. 
#' @return an \code{[n]} vector, containing binary predictors where the \code{0} group is sampled from
#' \eqn{Bern(balance/2)} and the \code{K-1} groups are sampled from \eqn{Bern(1 - balance/2)}.
#' @noRd
cb.sims.cate.simulate_covars_binary <- function(groups, balance=1) {
  covars <- numeric(length(groups))
  plow <- balance/2
  phigh <- 1 - plow
  
  for (i in 1:length(groups)) {
    if (groups[i] == 0) {
      covars[i] <- rbinom(1, 1, plow)
    } else {
      covars[i] <- rbinom(1, 1, phigh)
    }
  }
  
  return(covars)
}

#' Simulate Continuous Covariates
#' @importFrom stats rbinom rbeta
#' @param groups group assignments, as an \code{[n]} vector, which take values from \code{0} to \code{K} for \code{K} groups.
#' @param balance a parameter that governs the similarity between the continuous covariate distributions for the \code{0} group
#' against the other \code{K-1} groups.
#' @param alpha the alpha for sampling the \code{0} group, and the beta for sampling the other \code{K - 1} groups.
#' @param beta the beta for sampling the \code{0} group, and the alpha for sampling the other \code{K - 1} groups.
#' @param common a parameter which governs the shape of the common sampling distribution.
#' @return an \code{[n]} vector, containing binary predictors where the \code{0} group is sampled from
#' \eqn{Beta(alpha, beta)} with probability \code{1 - balance} and the \code{K-1} groups are sampled from \eqn{Beta(beta, alpha)} with probability \code{1 - balance},
#' and both groups are sampled from \eqn{Beta(common, common)} with probability \code{balance}..
#' @noRd
cb.sims.cate.simulate_covars <- function(groups, balance=1, alpha=2, beta=8, common=10) {
  coefs <- list()
  balance_id <- rbinom(length(groups), 1, balance)
  
  for (i in 1:length(groups)) {
    causal_cl <- groups[i]
    bal_id <- balance_id[i]
    
    if (bal_id == 1) {
      coefs[[i]] <- c(common, common)
    } else {
      if (causal_cl == 0) {
        coefs[[i]] <- c(alpha, beta)
      } else {
        coefs[[i]] <- c(beta, alpha)
      }
    }
  }
  
  covars <- numeric(length(groups))
  for (i in 1:length(coefs)) {
    alpha <- coefs[[i]][1]
    beta <- coefs[[i]][2]
    covars[i] <- 2 * rbeta(1, alpha, beta) - 1
  }
  
  return(covars)
}

#' Simulate Continuous Covariates for Multiple Classes
#' @importFrom stats rbinom rbeta
#' @param groups group assignments, as an \code{[n]} vector, which take values from \code{0} to \code{K-1} for \code{K} groups.
#' @param balance a parameter that governs the similarity between the continuous covariate distributions for the \code{0} group
#' against the other \code{K-1} groups.
#' @param alpha the alpha for sampling the \code{0} group, and the beta for sampling the other \code{K-1} groups.
#' @param beta the beta for sampling the \code{0} group, and the alpha for sampling the other \code{K-1} groups.
#' @param common a parameter which governs the shape of the common sampling distribution.
#' @return an \code{[n]} vector, containing continuous predictors where the \code{0} group is sampled from
#' \eqn{Beta(alpha, beta)} with probability \code{1 - balance} and the \code{K-1} groups are sampled from \eqn{Beta(beta, alpha)} with probability \code{1 - balance},
#' and all groups are sampled from \eqn{Beta(common, common)} with probability \code{balance}.
#' @noRd
cb.sims.cate.simulate_covars_multiclass <- function(groups, balance=1, alpha=2, beta=8, common=10) {
  coefs <- list()
  balance_id <- rbinom(length(groups), 1, balance)
  
  for (i in 1:length(groups)) {
    causal_cl <- groups[i]
    bal_id <- balance_id[i]
    
    if (bal_id == 1) {
      coefs[[i]] <- c(common, common)
    } else {
      if (causal_cl == 0) {
        coefs[[i]] <- c(alpha, beta)
      } else {
        coefs[[i]] <- c(beta, alpha)
      }
    }
  }
  
  covars <- numeric(length(groups))
  for (i in 1:length(coefs)) {
    alpha <- coefs[[i]][1]
    beta <- coefs[[i]][2]
    covars[i] <- 2 * rbeta(1, alpha, beta) - 1
  }
  
  return(covars)
}


#' Generate a Random Orthogonal Matrix
#'
#' Creates a random orthogonal matrix from the special orthogonal group SO(n),
#' which is the group of n×n orthogonal matrices with determinant 1.
#' This is equivalent to scipy's \code{special_ortho_group.rvs} function.
#'
#' @param p Integer specifying the dimension of the square matrix to generate.
#'
#' @return An p×p orthogonal matrix with determinant 1.
#'
#' @details
#' The function generates a random matrix with standard normal entries,
#' performs QR decomposition to obtain an orthogonal matrix,
#' and ensures the determinant is 1 by potentially flipping a column sign.
#'
#' @examples
#' # Generate a 3×3 random orthogonal matrix
#' R <- random_orthogonal_matrix(3)
#' 
#' # Verify orthogonality: R'R should be approximately the identity matrix
#' t(R) %*% R
#' 
#' # Verify determinant is 1
#' det(R)
#'
#' @importFrom stats rnorm
#' @noRd
cb.sims.random_rotation <- function(p) {
  # Generate a random matrix with standard normal entries
  M <- matrix(rnorm(p*p), nrow=p)
  
  # QR decomposition
  Q <- qr.Q(qr(M))
  
  # Ensure the determinant is 1 (special orthogonal group)
  if (det(Q) < 0) {
    # Flip the sign of the first column if determinant is negative
    Q[,1] <- -Q[,1]
  }
  
  return(Q)
}

#' Sigmoidal CATE Simulation
#' 
#' @param n the number of samples. Defaults to \code{100}.
#' @param p the number of dimensions. Defaults to \code{10}.
#' @param pi the balance between the classes, where samples will be from group 1
#' with probability \code{pi}, and group 2 with probability \code{1 - pi}. Defaults
#' to \code{0.5}.
#' @param balance a parameter governing the covariate similarity between the two groups. Defaults to \code{1}, which is the same covariate distributions for both groups.
#' @param eff_sz the conditional treatment effect size between the different groups, which governs the rotation in radians between the first and second group. Defaults to \code{1}.
#' @param covar_eff_sz A parameter which governs the covariate effect size with respect to the outcome. Defaults to \code{5}..
#' @param alpha the alpha for sampling the \code{0} group, and the beta for sampling the other \code{K - 1} groups.
#' @param beta the beta for sampling the \code{0} group, and the alpha for sampling the other \code{K - 1} groups.
#' @param common a parameter which governs the shape of the common sampling distribution.
#' @param pi the fraction of points which are sampled from a common distribution. Should be a number between \code{0} and \code{1}.
#' @param a the first parameter for the covariate/outcome relationship. Defaults to \code{2}.
#' @param b the second parameter for the covariate/outcome relationship. Defaults to \code{8}.
#' @param err the level of noise for the simulation. Defaults to \code{1}.
#' @param nbreaks the number of breakpoints for computing the expected outcome at a given covariate level
#' for each batch. Defaults to \code{200}.
#' @return a list, containing the following:
#' \item{Y}{an \code{[n, p]} matrix, containing the outcomes for each sample.}
#' \item{Ts}{an \code{[n, 1]} matrix, containing the group labels for each sample.}
#' \item{Xs}{an \code{[n, 1]} matrix, containing the covariate values for each sample.}
#' \item{Eps}{an \code{[n, 1]} matrix, containing the error for each sample.}
#' \item{Ytrue}{an \code{[nbreaks*2, 2]} matrix, containing the expected outcomes at a covariate level indicated by \code{Xtrue}.}
#' \item{Ttrue}{an \code{[nbreaks*2,1]} matrix, indicating the group/batch the expected outcomes and covariate breakpoints correspond to.}
#' \item{Xtrue}{an \code{[nbreaks*2, 1]} matrix, indicating the values of the covariate breakpoints for the theoretical expected outcome in \code{Ytrue}.}
#' \item{Group.Effect}{The group effect magnitude.}
#' \item{Covar.Effect}{The covariate effect magnitude.}
#' 
#' @references Eric W. Bridgeford, et al. "Learning Sources of Variability from High-Dimensional Observational Studies" arXiv (2025). 
#' 
#' @author Eric W. Bridgeford
#' 
#' @examples
#' 
#' library(causalBatch)
#' sim = cb.sims.cate.sigmoidal_sim()
#' 
#' @export
cb.sims.cate.sigmoidal_sim <- function(n=100, p=10, pi=0.5, balance=1, eff_sz=1, covar_eff_sz=5, alpha=2, beta=8, common=10, a=2, b=8, err=1, nbreaks=200) {
  rotation <- eff_sz * base::pi  # Angle of rotation of the second group
  rot_rescale <- cos(rotation)  # the rescaling factor for the rotation of the second group
  
  Ts <- rbinom(n, 1, pi)
  Xs <- cb.sims.cate.simulate_covars(Ts, balance=balance, alpha=alpha, beta=beta, common=common)
  y_base <- matrix(sigmoid(b*Xs), ncol=1)
  Bs <- a/((1:p)^1.5)
  Bs <- matrix(Bs, nrow=1)
  
  Ys_covar <- covar_eff_sz * (y_base %*% Bs)
  
  rot_vec <- rep(1, n)
  rot_vec[Ts == 0] <- rot_rescale
  R <- diag(rot_vec)
  
  Ys_covar <- R %*% (Ys_covar - covar_eff_sz/2 * matrix(rep(Bs, each=n), nrow=n))
  Ys_covar <- Ys_covar + covar_eff_sz/2 * matrix(rep(Bs, each=n), nrow=n)
  eps <- matrix(rnorm(n*p, sd=err), nrow=n, ncol=p)
  Ys <- Ys_covar + eps
  
  true_x <- seq(-1, 1, length.out=nbreaks)
  true_x <- c(true_x, true_x)
  true_y_base <- matrix(sigmoid(b*true_x), ncol=1)
  true_y_covar <- covar_eff_sz * (true_y_base %*% Bs)
  true_t <- c(rep(0, nbreaks), rep(1, nbreaks))
  rot_vec_true <- rep(1, 2*nbreaks)
  rot_vec_true[true_t == 0] <- rot_rescale
  R_true <- diag(rot_vec_true)
  true_y <- R_true %*% (true_y_covar - covar_eff_sz/2 * matrix(rep(Bs, each=2*nbreaks), nrow=2*nbreaks))
  true_y <- true_y + covar_eff_sz/2 * matrix(rep(Bs, each=2*nbreaks), nrow=2*nbreaks)
  
  out <- list(Ys=matrix(Ys, ncol=p), Ts=Ts, Xs=matrix(Xs, ncol=1), Eps=eps, Ytrue=true_y, 
              Ttrue=true_t, Xtrue=true_x, Group.Effect=eff_sz, Covar.Effect=covar_eff_sz)
  if (rotate) {
    # if desired, generate and apply a rotation matrix
    R <- cb.sims.random_rotation(p)
    out$Ys <- out$Ys %*% t(R)
    out$R <- R
  }
  return(out)
}

#' Non-monotone CATE Simulation
#' 
#' @param n the number of samples. Defaults to \code{100}.
#' @param p the number of dimensions. Defaults to \code{10}.
#' @param pi the balance between the classes, where samples will be from group 1
#' with probability \code{pi}, and group 2 with probability \code{1 - pi}. Defaults
#' to \code{0.5}.
#' @param balance a parameter governing the covariate similarity between the two groups. Defaults to \code{1}, which is the same covariate distributions for both groups.
#' @param eff_sz the conditional treatment effect size between the different groups, which governs the rotation in radians between the first and second group. Defaults to \code{1}.
#' @param covar_eff_sz A parameter which governs the covariate effect size with respect to the outcome. Defaults to \code{5}..
#' @param alpha the alpha for sampling the \code{0} group, and the beta for sampling the other \code{K - 1} groups.
#' @param beta the beta for sampling the \code{0} group, and the alpha for sampling the other \code{K - 1} groups.
#' @param common a parameter which governs the shape of the common sampling distribution.
#' @param pi the fraction of points which are sampled from a common distribution. Should be a number between \code{0} and \code{1}.
#' @param err the level of noise for the simulation. Defaults to \code{1}.
#' @param nbreaks the number of breakpoints for computing the expected outcome at a given covariate level
#' for each batch. Defaults to \code{200}.
#' @param rotate whether to apply a random rotation to the outcomes. Defaults to \code{TRUE}.
#' @return a list, containing the following:
#' \item{Y}{an \code{[n, p]} matrix, containing the outcomes for each sample.}
#' \item{Ts}{an \code{[n, 1]} matrix, containing the group labels for each sample.}
#' \item{Xs}{an \code{[n, 1]} matrix, containing the covariate values for each sample.}
#' \item{Eps}{an \code{[n, 1]} matrix, containing the error for each sample.}
#' \item{Ytrue}{an \code{[nbreaks*2, 2]} matrix, containing the expected outcomes at a covariate level indicated by \code{Xtrue}.}
#' \item{Ttrue}{an \code{[nbreaks*2,1]} matrix, indicating the group/batch the expected outcomes and covariate breakpoints correspond to.}
#' \item{Xtrue}{an \code{[nbreaks*2, 1]} matrix, indicating the values of the covariate breakpoints for the theoretical expected outcome in \code{Ytrue}.}
#' \item{Group.Effect}{The group effect magnitude.}
#' \item{R}{Optional argument returned if a rotation is requested for the rotation matrix applied.}
#' 
#' @references Eric W. Bridgeford, et al. "Learning Sources of Variability from High-Dimensional Observational Studies" arXiv (2025). 
#' 
#' @author Eric W. Bridgeford
#' 
#' @examples
#' 
#' library(causalBatch)
#' sim = cb.sims.cate.nonmonotone_sim()
#' 
#' @export
cb.sims.cate.nonmonotone_sim <- function(n=100, p=10, pi=0.5, balance=1, eff_sz=1, alpha=2,
                                         beta=8, common=10, err=1, nbreaks=200, rotate=TRUE) {
  Ts <- rbinom(n, 1, pi)
  Xs <- cb.sims.cate.simulate_covars(Ts, balance=balance, alpha=alpha, beta=beta, common=common)
  
  y_base <- matrix(Xs, ncol=1)
  Bs <- 2/((1:p)^1.5)
  Bs <- matrix(Bs, nrow=1)
  
  Ys_covar <- matrix(0, nrow=n, ncol=p)
  idx <- which(Xs >= -0.3 & Xs <= 0.3)
  Ys_covar[idx,] <- eff_sz * matrix(rep(Bs, length(idx)), nrow=length(idx), byrow=TRUE)
  idx_0 <- which(Ts == 0)
  Ys_covar[idx_0,] <- -Ys_covar[idx_0,]
  
  eps <- matrix(rnorm(n*p, sd=err), nrow=n, ncol=p)
  Ys <- Ys_covar + eps
  
  # true signal at a given x
  true_x <- seq(-1, 1, length.out=nbreaks)
  true_x <- c(true_x, true_x)
  true_y_covar <- matrix(0, nrow=2*nbreaks, ncol=p)
  idx <- which(true_x >= -0.3 & true_x <= 0.3)
  true_y_covar[idx,] <- eff_sz * matrix(rep(Bs, length(idx)), nrow=length(idx), byrow=TRUE)
  true_t <- c(rep(0, nbreaks), rep(1, nbreaks))
  idx_0 <- which(true_t == 0)
  true_y_covar[idx_0,] <- -true_y_covar[idx_0,]
  
  true_y <- true_y_covar
  
  out <- list(Ys=matrix(Ys, ncol=p), Ts=Ts, Xs=matrix(Xs, ncol=1), Eps=eps, 
              Ytrue=true_y, Ttrue=true_t, Xtrue=true_x, Group.Effect=eff_sz)
  if (rotate) {
    # if desired, generate and apply a rotation matrix
    R <- cb.sims.random_rotation(p)
    out$Ys <- out$Ys %*% t(R)
    out$R <- R
  }
  return(out)
}

#' K-class Sigmoidal CATE Simulation
#' 
#' @param n the number of samples. Defaults to \code{100}.
#' @param p the number of dimensions. Defaults to \code{10}.
#' @param pi the balance between the classes, where samples will be from group 1
#' with probability \code{pi}, and group 2 with probability \code{1 - pi}. Defaults
#' to \code{0.5}.
#' @param balance a parameter governing the covariate similarity between the two groups. Defaults to \code{1}, which is the same covariate distributions for both groups.
#' @param eff_sz the conditional treatment effect size between the different groups, which governs the rotation in radians between the first and second group. Defaults to \code{1}.
#' @param covar_eff_sz A parameter which governs the covariate effect size with respect to the outcome. Defaults to \code{5}..
#' @param alpha the alpha for sampling the \code{0} group, and the beta for sampling the other \code{K - 1} groups.
#' @param beta the beta for sampling the \code{0} group, and the alpha for sampling the other \code{K - 1} groups.
#' @param common a parameter which governs the shape of the common sampling distribution.
#' @param pi the fraction of points which are sampled from a common distribution. Should be a number between \code{0} and \code{1}.
#' @param a the first parameter for the covariate/outcome relationship. Defaults to \code{2}.
#' @param b the second parameter for the covariate/outcome relationship. Defaults to \code{8}.
#' @param err the level of noise for the simulation. Defaults to \code{1}.
#' @param nbreaks the number of breakpoints for computing the expected outcome at a given covariate level
#' for each batch. Defaults to \code{200}.
#' @param K the number of classes. Defaults to \code{3}.
#' @param rotate whether to apply a random rotation to the outcomes. Defaults to \code{TRUE}.
#' @return a list, containing the following:
#' \item{Y}{an \code{[n, p]} matrix, containing the outcomes for each sample.}
#' \item{Ts}{an \code{[n, 1]} matrix, containing the group labels for each sample.}
#' \item{Xs}{an \code{[n, 1]} matrix, containing the covariate values for each sample.}
#' \item{Eps}{an \code{[n, 1]} matrix, containing the error for each sample.}
#' \item{Ytrue}{an \code{[nbreaks*2, 2]} matrix, containing the expected outcomes at a covariate level indicated by \code{Xtrue}.}
#' \item{Ttrue}{an \code{[nbreaks*2,1]} matrix, indicating the group/batch the expected outcomes and covariate breakpoints correspond to.}
#' \item{Xtrue}{an \code{[nbreaks*2, 1]} matrix, indicating the values of the covariate breakpoints for the theoretical expected outcome in \code{Ytrue}.}
#' \item{Group.Effect}{The group effect magnitude.}
#' \item{Covar.Effect}{The covariate effect magnitude.}
#' \item{K}{The total number of classes.}
#' \item{R}{Optional argument returned if a rotation is requested for the rotation matrix applied.}
#' 
#' @references Eric W. Bridgeford, et al. "Learning Sources of Variability from High-Dimensional Observational Studies" arXiv (2025). 
#' 
#' @author Eric W. Bridgeford
#' 
#' @examples
#' 
#' library(causalBatch)
#' sim = cb.sims.cate.kclass_rotation_sim()
#' 
#' @export
cb.sims.cate.kclass_sigmoidal_sim <- function(n=100, p=10, pi=0.5, balance=1, eff_sz=1, covar_eff_sz=5, 
                                             alpha=2, beta=8, common=10, a=2, b=8, err=1, nbreaks=200, K=3,
                                             rotate=TRUE) {
  rotation <- eff_sz * base::pi  # Angle of rotation of the second group
  rot_rescale <- cos(rotation)  # the rescaling factor for the rotation of the second group
  
  probs <- c(pi, rep((1-pi)/(K-1), K-1))
  Ts <- sample(0:(K-1), n, replace=TRUE, prob=probs)
  Xs <- cb.sims.cate.simulate_covars_multiclass(Ts, balance=balance, alpha=alpha, beta=beta, common=common)
  y_base <- matrix(sigmoid(b*Xs), ncol=1)
  Bs <- a/((1:p)^1.1)
  Bs <- matrix(Bs, nrow=1)
  
  Ys_covar <- covar_eff_sz * (y_base %*% Bs)
  
  rot_vec <- rep(1, n)
  rot_vec[Ts == 0] <- rot_rescale
  R <- diag(rot_vec)
  
  Ys_covar <- R %*% (Ys_covar - covar_eff_sz/2 * matrix(rep(Bs, each=n), nrow=n))
  Ys_covar <- Ys_covar + covar_eff_sz/2 * matrix(rep(Bs, each=n), nrow=n)
  eps <- matrix(rnorm(n*p, sd=err), nrow=n, ncol=p)
  Ys <- Ys_covar + eps
  
  # true signal at a given x
  Ntrue <- K*nbreaks
  true_x <- seq(-1, 1, length.out=nbreaks)
  true_x <- rep(true_x, K)
  true_y_base <- matrix(sigmoid(b*true_x), ncol=1)
  true_y_covar <- covar_eff_sz * (true_y_base %*% Bs)
  true_t <- as.vector(sapply(0:(K-1), function(k) rep(k, nbreaks)))
  rot_vec_true <- rep(1, Ntrue)
  rot_vec_true[true_t == 0] <- rot_rescale
  R_true <- diag(rot_vec_true)
  true_y <- R_true %*% (true_y_covar - covar_eff_sz/2 * matrix(rep(Bs, each=Ntrue), nrow=Ntrue))
  true_y <- true_y + covar_eff_sz/2 * matrix(rep(Bs, each=Ntrue), nrow=Ntrue)
  
  out <- list(Ys=matrix(Ys, ncol=p), Ts=Ts, Xs=matrix(Xs, ncol=1), Eps=eps, Ytrue=true_y, Ttrue=true_t, Xtrue=true_x, 
              Group.Effect=eff_sz, Covar.Effect=covar_eff_sz, K=K)
  if (rotate) {
    # if desired, generate and apply a rotation matrix
    R <- cb.sims.random_rotation(p)
    out$Ys <- out$Ys %*% t(R)
    out$R <- R
  }
  return(out)
}

#' Simulate Data with Heteroskedastic Conditional Average Treatment Effects
#' @importFrom stats rbinom rnorm
#' @param n number of samples to generate
#' @param p number of dimensions for the response variable
#' @param pi probability of assignment to treatment group 1
#' @param balance a parameter that governs the similarity between the covariate distributions
#' @param eff_sz effect size parameter controlling the heteroskedasticity between groups
#' @param covar_eff_sz effect size parameter for the covariate influence
#' @param alpha the alpha parameter for beta distribution of control group
#' @param beta the beta parameter for beta distribution of control group
#' @param common parameter governing the shape of the common sampling distribution
#' @param a scaling factor for response variable generation
#' @param b scaling factor for sigmoid transformation
#' @param err standard deviation for the error term
#' @param nbreaks number of points to use for generating the true signal
#' @param rotate whether to apply a random rotation to the outcomes. Defaults to \code{TRUE}.
#' @return A list containing:
#'   \item{Ys}{response matrix of size [n, p]}
#'   \item{Ts}{treatment assignment vector of size [n]}
#'   \item{Xs}{covariate vector of size [n]}
#'   \item{Eps}{error matrix of size [n, p]}
#'   \item{Ytrue}{true response matrix for evaluation}
#'   \item{Ttrue}{true treatment vector for evaluation}
#'   \item{Xtrue}{true covariate vector for evaluation}
#'   \item{Group.Effect}{the effect size parameter}
#'   \item{Covar.Effect}{the covariate effect size parameter}
#' \item{R}{Optional argument returned if a rotation is requested for the rotation matrix applied.}
#' 
#' @references Eric W. Bridgeford, et al. "Learning Sources of Variability from High-Dimensional Observational Studies" arXiv (2025). 
#' 
#' @author Eric W. Bridgeford
#' 
#' @examples
#' 
#' library(causalBatch)
#' sim = cb.sims.cate.heteroskedastic_sim()
#' 
#' @export
cb.sims.cate.heteroskedastic_sim <- function(n=100, p=10, pi=0.5, balance=1, eff_sz=1, covar_eff_sz=3, 
                                             alpha=2, beta=8, common=10, a=2, b=8, err=0.5, nbreaks=200,
                                             rotate = TRUE) {
  # First get base data from sigmoidal_sim_cate with zero effect size
  result <- cb.sims.cate.sigmoidal_sim(n=n, p=p, pi=pi, balance=balance, eff_sz=0, 
                                       covar_eff_sz=covar_eff_sz, alpha=alpha, beta=beta, 
                                       common=common, a=a, b=b, err=err, nbreaks=nbreaks, rotate=FALSE)
  
  Ys <- result$Ys
  Ts <- result$Ts
  
  # Add heteroskedastic noise to control group
  idx <- which(Ts == 0)
  hetero_noise <- matrix(rnorm(length(idx)*p, sd=sqrt(2*eff_sz)), nrow=length(idx), ncol=p)
  Ys[idx,] <- Ys[idx,] + hetero_noise
  
  out <- list(Ys=matrix(Ys, ncol=p), 
              Ts=Ts, 
              Xs=matrix(result$Xs, ncol=1), 
              Eps=result$Eps, 
              Ytrue=result$Ytrue, 
              Ttrue=result$Ttrue, 
              Xtrue=matrix(result$Xtrue, ncol=1), 
              Group.Effect=eff_sz, 
              Covar.Effect=covar_eff_sz)
  if (rotate) {
    # if desired, generate and apply a rotation matrix
    R <- cb.sims.random_rotation(p)
    out$Ys <- out$Ys %*% t(R)
    out$R <- R
  }
  return(out)
}