#' Adjust for batch effects using an empirical Bayes framework
#'
#' ComBat allows users to adjust for batch effects in datasets where the batch covariate is known, using methodology
#' described in Johnson et al. 2007. It uses either parametric or non-parametric empirical Bayes frameworks for adjusting data for
#' batch effects.  Users are returned an expression matrix that has been corrected for batch effects. The input
#' data are assumed to be cleaned and normalized before batch effect removal.
#'
#' @param dat Genomic measure matrix (dimensions probe x sample) - for example, expression matrix
#' @param batch {Batch covariate (only one batch allowed)}
#' @param mod Model matrix for outcome of interest and other covariates besides batch
#' @param par.prior (Optional) TRUE indicates parametric adjustments will be used, FALSE indicates non-parametric adjustments will be used
#' @param prior.plots (Optional) TRUE give prior plots with black as a kernel estimate of the empirical batch effect density and red as the parametric
#' @param mean.only (Optional) FALSE If TRUE ComBat only corrects the mean of the batch effect (no scale adjustment)
#' @param ref.batch (Optional) NULL If given, will use the selected batch as a reference for batch adjustment.
#' @param BPPARAM (Optional) BiocParallelParam for parallel operation
#'
#' @return data A probe x sample genomic measure matrix, adjusted for batch effects.
#'
#' @importFrom graphics lines par
#' @importFrom stats cor density dnorm model.matrix pf ppoints prcomp predict
#' qgamma qnorm qqline qqnorm qqplot smooth.spline var
#' @importFrom utils read.delim
#'
#' @examples
#' library(bladderbatch)
#' data(bladderdata)
#' dat <- bladderEset[1:50,]
#'
#' pheno = pData(dat)
#' edata = exprs(dat)
#' batch = pheno$batch
#' mod = model.matrix(~as.factor(cancer), data=pheno)
#'
#' # parametric adjustment
#' combat_edata1 = ComBat(dat=edata, batch=batch, mod=NULL, par.prior=TRUE, prior.plots=FALSE)
#'
#' # non-parametric adjustment, mean-only version
#' combat_edata2 = ComBat(dat=edata, batch=batch, mod=NULL, par.prior=FALSE, mean.only=TRUE)
#'
#' # reference-batch version, with covariates
#' combat_edata3 = ComBat(dat=edata, batch=batch, mod=mod, par.prior=TRUE, ref.batch=3)
#'
#'

ComBat_return_model <- function (dat, batch, mod = NULL, par.prior = TRUE, prior.plots = FALSE, 
          mean.only = FALSE, ref.batch = NULL, BPPARAM = bpparam("SerialParam")) 
{
  if (length(dim(batch)) > 1) {
    stop("This version of ComBat only allows one batch variable")
  }
  dat <- as.matrix(dat)
  batch <- as.factor(batch)
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
  a.prior <- apply(delta.hat, 1, sva:::aprior)
  b.prior <- apply(delta.hat, 1, sva:::bprior)
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
        temp <- sva:::it.sol(s.data[, batches[[i]]], gamma.hat[i, 
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
  return(list(Corrected=bayesdata, Bayesdata=s.data, Var=var.pooled, Grand.mean=grand.mean,
              B.hat=B.hat, Design=batch.design, Gamma=gamma.star, Delta=delta.star, N=n.batches, 
              Stand=stand.mean))
}

apply_batch_adjustment <- function(data, batches, design, gamma, delta, N, var.pooled, grand.mean, B.hat, mod) {
  n.array <- dim(data)[2]
  n.batch <- length(unique(batches))
  
  design <- as.matrix(design)
  stand.mean <- t(grand.mean) %*% t(rep(1, n.array))

  design.mod <- cbind(design, mod)
  design.mod <- design.mod[,which(!sapply(1:dim(design.mod)[2], function(j) {all(design.mod[,j] == 1)}))]

  if (!is.null(design.mod)) {
    tmp <- as.matrix(design.mod)
    tmp[, c(1:n.batch)] <- 0
    stand.mean <- stand.mean + t(tmp %*% B.hat)
  }
  
  s.data <- (data - stand.mean)/(sqrt(var.pooled) %*% t(rep(1, 
                                                           n.array)))
  j <- 1
  bayesdata <- s.data
  for (i in batches) {
    bayesdata[, i] <- (bayesdata[, i] - t(design[i,] %*% gamma))/(sqrt(delta[j, ]) %*% t(rep(1, length(i))))
    j <- j + 1
  }
  
  bayesdata <- (bayesdata * (sqrt(var.pooled) %*% t(rep(1, n.array)))) + stand.mean
  
  return(t(bayesdata))
}