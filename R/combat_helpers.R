#' Adjust for batch effects using an empirical Bayes framework
#'
#' ComBat allows users to adjust for batch effects in datasets where the batch covariate is known, using methodology
#' described in Johnson et al. 2007. It uses either parametric or non-parametric empirical Bayes frameworks for adjusting data for
#' batch effects.  Users are returned an expression matrix that has been corrected for batch effects. The input
#' data are assumed to be cleaned and normalized before batch effect removal.
#' 
#' Note: this code is adapted directly from the \code{\link[sva]{ComBat}} algorithm featured in the `sva` package, and is not intended for standalone use.
#'
#' @param dat Genomic measure matrix (sample x dimensions probe) - for example, expression matrix
#' @param batch {Batch covariate (only one batch allowed)}
#' @param mod Model matrix for outcome of interest and other covariates besides batch
#' @param par.prior (Optional) TRUE indicates parametric adjustments will be used, FALSE indicates non-parametric adjustments will be used
#' @param prior.plots (Optional) TRUE give prior plots with black as a kernel estimate of the empirical batch effect density and red as the parametric
#' @param mean.only (Optional) FALSE If TRUE ComBat only corrects the mean of the batch effect (no scale adjustment)
#' @param ref.batch (Optional) NULL If given, will use the selected batch as a reference for batch adjustment.
#' @param BPPARAM (Optional) BiocParallelParam for parallel operation
#' 
#' @return a list containing:
#' \itemize{
#'    \item{\code{Ys.corrected}} batch effect corrected data.
#'    \item{\code{Model}} the learned batch effect correction model.
#' }
#'
#' @importFrom graphics lines par
#' @importFrom stats cor density dnorm model.matrix pf ppoints prcomp predict
#' qgamma qnorm qqline qqnorm qqplot smooth.spline var
#' @importFrom utils read.delim
#' @importFrom genefilter rowVars
#' @import sva
#' @importFrom BiocParallel bplapply bpparam
#' @references W Evan Johnson, et al. "Adjusting batch effects in microarray expression data using empirical Bayes methods" Biostatistics (2007).
#' @references Leek JT, Johnson WE, Parker HS, Fertig EJ, Jaffe AE, Zhang Y, Storey JD, Torres LC (2024). sva: Surrogate Variable Analysis. R package version 3.52.0.
#'
#' @noRd
cb.learn.fit_cComBat <- function (dat, batch, mod = NULL, par.prior = TRUE, prior.plots = FALSE, 
                                 mean.only = FALSE, ref.batch = NULL, BPPARAM = bpparam("SerialParam")) 
{
  if (length(dim(batch)) > 1) {
    stop("This version of ComBat only allows one batch variable")
  }
  dat <- as.matrix(t(dat))
  batch <- as.factor(batch)
  batch.levels <- levels(batch)
  zero.rows.lst <- lapply(levels(batch), function(batch_level) {
    if (sum(batch == batch_level) > 1) {
      return(which(apply(dat[, batch == batch_level], 1, 
                         function(x) {
                           var(x) == 0
                         })))
    }
    else {
      return(which(rep(1, 3) == 2))
    }
  })
  zero.rows <- Reduce(union, zero.rows.lst)
  keep.rows <- setdiff(1:nrow(dat), zero.rows)
  if (length(zero.rows) > 0) {
    cat(sprintf("Found %d genes with uniform expression within a single batch (all zeros); these will not be adjusted for batch.\n", 
                length(zero.rows)))
    dat.orig <- dat
    dat <- dat[keep.rows, ]
  }
  if (any(table(batch) == 1)) {
    mean.only = TRUE
  }
  if (mean.only == TRUE) {
    message("Using the 'mean only' version of ComBat")
  }
  batchmod <- model.matrix(~-1 + batch)
  if (!is.null(ref.batch)) {
    if (!(ref.batch %in% levels(batch))) {
      stop("reference level ref.batch is not one of the levels of the batch variable")
    }
    message("Using batch =", ref.batch, "as a reference batch (this batch won't change)")
    ref <- which(levels(as.factor(batch)) == ref.batch)
    batchmod[, ref] <- 1
  } else {
    ref <- NULL
  }
  
  n.batch <- nlevels(batch)
  batches <- list()
  for (i in 1:n.batch) {
    batches[[i]] <- which(batch == levels(batch)[i])
  }
  n.batches <- sapply(batches, length)
  if (any(n.batches == 1)) {
    mean.only = TRUE
    message("Note: one batch has only one sample, setting mean.only=TRUE")
  }
  n.array <- sum(n.batches)
  design <- cbind(batchmod, mod)
  check <- apply(design, 2, function(x) all(x == 1))
  if (!is.null(ref)) {
    check[ref] <- FALSE
  }
  design <- as.matrix(design[, !check])
  message("Adjusting for", ncol(design) - ncol(batchmod), "covariate(s) or covariate level(s)")
  if (qr(design)$rank < ncol(design)) {
    if (ncol(design) == (n.batch + 1)) {
      stop("The covariate is confounded with batch! Remove the covariate and rerun ComBat")
    }
    if (ncol(design) > (n.batch + 1)) {
      if ((qr(design[, -c(1:n.batch)])$rank < ncol(design[, 
                                                          -c(1:n.batch)]))) {
        stop("The covariates are confounded! Please remove one or more of the covariates so the design is not confounded")
      } else {
        stop("At least one covariate is confounded with batch! Please remove confounded covariates and rerun ComBat")
      }
    }
  }
  NAs <- any(is.na(dat))
  if (NAs) {
    message(c("Found", sum(is.na(dat)), "Missing Data Values"), 
            sep = " ")
  }
  
  if (!NAs) {
    B.hat <- solve(crossprod(design), tcrossprod(t(design), 
                                                 as.matrix(dat)))
  } else {
    B.hat <- apply(dat, 1, Beta.NA, design)
  }
  if (!is.null(ref.batch)) {
    grand.mean <- t(B.hat[ref, ])
  } else {
    grand.mean <- crossprod(n.batches/n.array, B.hat[1:n.batch, 
    ])
  }
  if (!NAs) {
    if (!is.null(ref.batch)) {
      ref.dat <- dat[, batches[[ref]]]
      var.pooled <- ((ref.dat - t(design[batches[[ref]], 
      ] %*% B.hat))^2) %*% rep(1/n.batches[ref], n.batches[ref])
    } else {
      var.pooled <- ((dat - t(design %*% B.hat))^2) %*% 
        rep(1/n.array, n.array)
    }
  } else {
    if (!is.null(ref.batch)) {
      ref.dat <- dat[, batches[[ref]]]
      var.pooled <- rowVars(ref.dat - t(design[batches[[ref]], 
      ] %*% B.hat), na.rm = TRUE)
    } else {
      var.pooled <- rowVars(dat - t(design %*% B.hat), 
                            na.rm = TRUE)
    }
  }
  stand.mean <- t(grand.mean) %*% t(rep(1, n.array))
  if (!is.null(design)) {
    tmp <- design
    tmp[, c(1:n.batch)] <- 0
    stand.mean <- stand.mean + t(tmp %*% B.hat)
  }
  s.data <- (dat - stand.mean)/(sqrt(var.pooled) %*% t(rep(1, 
                                                           n.array)))
  
  batch.design <- design[, 1:n.batch]
  if (!NAs) {
    gamma.hat <- solve(crossprod(batch.design), tcrossprod(t(batch.design), 
                                                           as.matrix(s.data)))
  } else {
    gamma.hat <- apply(s.data, 1, Beta.NA, batch.design)
  }
  delta.hat <- NULL
  for (i in batches) {
    if (mean.only == TRUE) {
      delta.hat <- rbind(delta.hat, rep(1, nrow(s.data)))
    } else {
      delta.hat <- rbind(delta.hat, rowVars(s.data[, i], 
                                            na.rm = TRUE))
    }
  }
  gamma.bar <- rowMeans(gamma.hat)
  t2 <- rowVars(gamma.hat)
  a.prior <- apply(delta.hat, 1, aprior)
  b.prior <- apply(delta.hat, 1, bprior)
  if (prior.plots && par.prior) {
    old_pars <- par(no.readonly = TRUE)
    on.exit(par(old_pars))
    par(mfrow = c(2, 2))
    tmp <- density(gamma.hat[1, ])
    plot(tmp, type = "l", main = expression(paste("Density Plot of First Batch ", 
                                                  hat(gamma))))
    xx <- seq(min(tmp$x), max(tmp$x), length = 100)
    lines(xx, dnorm(xx, gamma.bar[1], sqrt(t2[1])), col = 2)
    qqnorm(gamma.hat[1, ], main = expression(paste("Normal Q-Q Plot of First Batch ", 
                                                   hat(gamma))))
    qqline(gamma.hat[1, ], col = 2)
    tmp <- density(delta.hat[1, ])
    xx <- seq(min(tmp$x), max(tmp$x), length = 100)
    tmp1 <- list(x = xx, y = dinvgamma(xx, a.prior[1], b.prior[1]))
    plot(tmp, typ = "l", ylim = c(0, max(tmp$y, tmp1$y)), 
         main = expression(paste("Density Plot of First Batch ", 
                                 hat(delta))))
    lines(tmp1, col = 2)
    invgam <- 1/qgamma(1 - ppoints(ncol(delta.hat)), a.prior[1], 
                       b.prior[1])
    qqplot(invgam, delta.hat[1, ], main = expression(paste("Inverse Gamma Q-Q Plot of First Batch ", 
                                                           hat(delta))), ylab = "Sample Quantiles", xlab = "Theoretical Quantiles")
    lines(c(0, max(invgam)), c(0, max(invgam)), col = 2)
  }
  gamma.star <- delta.star <- matrix(NA, nrow = n.batch, ncol = nrow(s.data))
  if (par.prior) {
    
    results <- bplapply(1:n.batch, function(i) {
      if (mean.only) {
        gamma.star <- postmean(gamma.hat[i, ], gamma.bar[i], 
                               1, 1, t2[i])
        delta.star <- rep(1, nrow(s.data))
      }
      else {
        temp <- it.sol(s.data[, batches[[i]]], gamma.hat[i, 
        ], delta.hat[i, ], gamma.bar[i], t2[i], a.prior[i], 
        b.prior[i])
        gamma.star <- temp[1, ]
        delta.star <- temp[2, ]
      }
      list(gamma.star = gamma.star, delta.star = delta.star)
    }, BPPARAM = BPPARAM)
    for (i in 1:n.batch) {
      gamma.star[i, ] <- results[[i]]$gamma.star
      delta.star[i, ] <- results[[i]]$delta.star
    }
  } else {
    results <- bplapply(1:n.batch, function(i) {
      if (mean.only) {
        delta.hat[i, ] = 1
      }
      temp <- int.eprior(as.matrix(s.data[, batches[[i]]]), 
                         gamma.hat[i, ], delta.hat[i, ])
      list(gamma.star = temp[1, ], delta.star = temp[2, 
      ])
    }, BPPARAM = BPPARAM)
    for (i in 1:n.batch) {
      gamma.star[i, ] <- results[[i]]$gamma.star
      delta.star[i, ] <- results[[i]]$delta.star
    }
  }
  if (!is.null(ref.batch)) {
    gamma.star[ref, ] <- 0
    delta.star[ref, ] <- 1
  }
  
  bayesdata <- s.data
  j <- 1
  for (i in batches) {
    bayesdata[, i] <- (bayesdata[, i] - t(batch.design[i, 
    ] %*% gamma.star))/(sqrt(delta.star[j, ]) %*% t(rep(1, 
                                                        n.batches[j])))
    j <- j + 1
  }
  bayesdata <- (bayesdata * (sqrt(var.pooled) %*% t(rep(1, 
                                                        n.array)))) + stand.mean
  
  if (!is.null(ref.batch)) {
    bayesdata[, batches[[ref]]] <- dat[, batches[[ref]]]
  }
  if (length(zero.rows) > 0) {
    dat.orig[keep.rows, ] <- bayesdata
    bayesdata <- dat.orig
  }
  return(list(Ys.corrected=t(bayesdata), Model=list(Var=var.pooled, Grand.mean=grand.mean,
                                                    B.hat=B.hat, Gamma=gamma.star, Delta=delta.star,
                                                    Levels=batch.levels)))
}

#' Density of inverse gamma distribution
#' 
#' Direct import of code from \code{\link[sva]{sva}} package.
#' @references Leek JT, Johnson WE, Parker HS, Fertig EJ, Jaffe AE, Zhang Y, Storey JD, Torres LC (2024). sva: Surrogate Variable Analysis. R package version 3.52.0.
#' @noRd
dinvgamma <- utils::getFromNamespace("dinvgamma", "sva")

#' Monte carlo integration
#' 
#' Direct import of code from \code{\link[sva]{sva}} package.
#' @references Leek JT, Johnson WE, Parker HS, Fertig EJ, Jaffe AE, Zhang Y, Storey JD, Torres LC (2024). sva: Surrogate Variable Analysis. R package version 3.52.0.
#' @noRd
int.eprior <- utils::getFromNamespace("int.eprior", "sva")

#' Fit LS models in presence of missing values
#' 
#' Direct import of code from \code{\link[sva]{sva}} package.
#' @references Leek JT, Johnson WE, Parker HS, Fertig EJ, Jaffe AE, Zhang Y, Storey JD, Torres LC (2024). sva: Surrogate Variable Analysis. R package version 3.52.0.
#' @noRd
Beta.NA <- utils::getFromNamespace("Beta.NA", "sva")

#' Postmean function for hyper-priors
#' 
#' Direct import of code from \code{\link[sva]{sva}} package.
#' @references Leek JT, Johnson WE, Parker HS, Fertig EJ, Jaffe AE, Zhang Y, Storey JD, Torres LC (2024). sva: Surrogate Variable Analysis. R package version 3.52.0.
#' @noRd
postmean <- utils::getFromNamespace("postmean", "sva")

#' Aprior function for hyper-priors
#' 
#' Direct import of code from \code{\link[sva]{sva}} package.
#' @references Leek JT, Johnson WE, Parker HS, Fertig EJ, Jaffe AE, Zhang Y, Storey JD, Torres LC (2024). sva: Surrogate Variable Analysis. R package version 3.52.0.
#' @noRd
aprior <- utils::getFromNamespace("aprior", "sva")

#' Bprior function for hyper-priors
#' 
#' Direct import of code from \code{\link[sva]{sva}} package.
#' @references Leek JT, Johnson WE, Parker HS, Fertig EJ, Jaffe AE, Zhang Y, Storey JD, Torres LC (2024). sva: Surrogate Variable Analysis. R package version 3.52.0.
#' @noRd
bprior <- utils::getFromNamespace("bprior", "sva")

#' Uses EM to find batch adjustments
#' 
#' Direct import of code from \code{\link[sva]{sva}} package.
#' @references Leek JT, Johnson WE, Parker HS, Fertig EJ, Jaffe AE, Zhang Y, Storey JD, Torres LC (2024). sva: Surrogate Variable Analysis. R package version 3.52.0.
#' @noRd
it.sol <- utils::getFromNamespace("it.sol", "sva")

#' Adjust for batch effects using an empirical Bayes framework
#'
#' ComBat allows users to adjust for batch effects in datasets where the batch covariate is known, using methodology
#' described in Johnson et al. 2007. It uses either parametric or non-parametric empirical Bayes frameworks for adjusting data for
#' batch effects.  Users are returned an expression matrix that has been corrected for batch effects. The input
#' data are assumed to be cleaned and normalized before batch effect removal.
#' 
#' Note: this code is adapted directly from the \code{\link[sva]{ComBat}} algorithm featured in the `sva` package.
#'
#' @param Ys an \code{[n, d]} matrix, for the outcome variables with \code{n} samples in \code{d} dimensions.
#' @param Ts \code{[n]} the labels of the samples, with \code{K < n} levels, as a factor variable.
#' @param Xs \code{[n, r]} the \code{r} covariates/confounding variables, for each of the \code{n} samples, as a data frame with named columns.
#' @param Model a list containing the following parameters:
#' \itemize{
#'    \item{\code{Var}} the pooled variance
#'    \item{\code{Grand.mean}} the overall mean of the data
#'    \item{\code{B.hat}} the fit regression coefficients
#'    \item{\code{Gamma}} additive batch effects
#'    \item{\code{Delta}} multiplicative batch effects
#'    \item{\code{Levels}} the order of levels for each batch
#'    \item{\code{Covar.Mod}} the covariate model for adjustment
#' }
#' This model is output after fitting with \code{\link{cb.correct.matching_cComBat}}.
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
cb.correct.apply_cComBat <- function(Ys, Ts, Xs, Model) {
  Ys <- t(Ys)
  n.array <- dim(Ys)[2]
  n.batch <- length(Model$Levels)
  
  batches <- lapply(Model$Levels, function(batch) {
    which(Ts == batch)
  })
  
  for (batch in Ts) {
    if (!(batch %in% Model$Levels)) {
      stop("You have out-of-sample data from batches not in the fit model.")
    }
  }
  
  design.batch <- as.matrix(ohe(Ts, levels=Model$Levels)$ohe)
  stand.mean <- t(Model$Grand.mean) %*% t(rep(1, n.array))
  
  if (!is.null(Model$Covar.Mod)) {
    mod.mtx <- model.matrix(as.formula(sprintf("~%s", Model$Covar.Mod)), data=Xs)
  } else {
    mod.mtx <- NULL
  }
  design.mod <- cbind(design.batch, mod.mtx)
  design.mod <- design.mod[,which(!sapply(1:dim(design.mod)[2], function(j) {all(design.mod[,j] == 1)}))]
  
  if (!is.null(design.mod)) {
    tmp <- as.matrix(design.mod)
    tmp[, c(1:n.batch)] <- 0
    stand.mean <- stand.mean + t(tmp %*% Model$B.hat)
  }
  
  s.data <- (Ys - stand.mean)/(sqrt(Model$Var) %*% t(rep(1, n.array)))
  j <- 1
  bayesdata <- s.data
  for (i in batches) {
    bayesdata[, i] <- (bayesdata[, i] - t(design.batch[i,] %*% Model$Gamma))/(sqrt(Model$Delta[j, ]) %*% t(rep(1, length(i))))
    j <- j + 1
  }
  
  bayesdata <- (bayesdata * (sqrt(Model$Var) %*% t(rep(1, n.array)))) + stand.mean
  
  return(t(bayesdata))
}