is.cran <- testthat:::on_cran()
ncores <- ifelse(is.cran, 1, parallel::detectCores())

# Setup parameters for test cases
R <- ifelse(is.cran, 100, 200)
nreps.null <- ifelse(is.cran, 50, 100)
alpha <- ifelse(is.cran, 0.25, 0.15)
eps <- ifelse(is.cran, 0.2, 0.1)
eff_sz <- 3
n <- ifelse(is.cran, 100, 150)


test_that("reject null when balanced and ATE", {
  sim.high <- cb.sims.sim_linear(n=n, eff_sz=eff_sz)
  
  res <- cb.detect.caus_cdcorr(sim.high$Ys, sim.high$Ts, sim.high$Xs,
                               num.threads = ncores,
                               R=R)
  expect_true(res$Test$p.value < alpha)
})

test_that("reject null when imbalanced and ATE", {
  sim.low <- cb.sims.sim_linear(n=n, eff_sz=eff_sz, unbalancedness=1.5)
  
  res <- cb.detect.caus_cdcorr(sim.low$Ys, sim.low$Ts, sim.low$Xs,
                               num.threads=ncores, R=R)
  expect_true(res$Test$p.value < alpha)
})

test_that("typically fails to reject null when balanced and ATE=0", {
  skip_on_cran()
  pwr_test <- sapply(1:nreps.null, function(i) {
    sim.high <- cb.sims.sim_linear(n=n, eff_sz=0)
    
    res <- cb.detect.caus_cdcorr(sim.high$Ys, sim.high$Ts, sim.high$Xs,
                                 num.threads = ncores, R=R)
    res$Test$p.value < alpha - eps
  })
  
  # should falsely reject null in favor of alternative at a rate ~alpha + small 
  # margin of error
  expect_true(mean(pwr_test) < alpha + eps)
})

test_that("typically fails to reject null when imbalanced and ATE=0", {
  skip_on_cran()
  pwr_test <- sapply(1:nreps.null, function(i) {
    sim.low <- cb.sims.sim_linear(n=n, eff_sz=0, unbalancedness = 1.5)
    
    res <- cb.detect.caus_cdcorr(sim.low$Ys, sim.low$Ts, sim.low$Xs,
                                 num.threads = ncores, R=R)
    res$Test$p.value < alpha - eps
  })
  
  expect_true(mean(pwr_test) < alpha + eps)
})

test_that("works when covariates are given as a vector", {
  sim.high <- cb.sims.sim_linear(n=n, eff_sz=eff_sz)
  
  res <- cb.detect.caus_cdcorr(sim.high$Ys, sim.high$Ts, as.vector(sim.high$Xs), R=R)
  expect_true(res$Test$p.value < alpha + eps)
})

test_that("works with >1 covariate levels", {
  sim.high <- cb.sims.sim_linear(n=n, unbalancedness=1, eff_sz=eff_sz)
  sim.high$Xs <- cbind(sim.high$Xs, runif(n))
  
  res <- cb.detect.caus_cdcorr(sim.high$Ys, sim.high$Ts, sim.high$Xs,
                               num.threads = ncores, R=R)
  expect_true(res$Test$p.value < alpha + eps)
})

test_that("rejects null with non-linearities and ATE", {
  sim.high <- cb.sims.sim_sigmoid(n=n, eff_sz=eff_sz)
  
  res <- cb.detect.caus_cdcorr(sim.high$Ys, sim.high$Ts, sim.high$Xs,
                               num.threads = ncores, R=R)
  expect_true(res$Test$p.value < alpha + eps)
  
  sim.low <- cb.sims.sim_sigmoid(n=n, eff_sz=eff_sz, unbalancedness=1.5)
  
  res <- cb.detect.caus_cdcorr(sim.low$Ys, sim.low$Ts, sim.low$Xs,
                               num.threads = ncores, R=R)
  expect_true(res$Test$p.value < alpha + eps)
})

test_that("fails to reject null with non-linearities and no ATE", {
  skip_on_cran()
  # balanced
  pwr_test.bal <- sapply(1:nreps.null, function(i) {
    sim.high <- cb.sims.sim_sigmoid(n=n, eff_sz=0)
    
    res <- cb.detect.caus_cdcorr(sim.high$Ys, sim.high$Ts, sim.high$Xs,
                                 num.threads = ncores, R=R)
    res$Test$p.value < alpha - eps
  })
  expect_true(mean(pwr_test.bal) < alpha + eps)
  
  # imbalanced
  pwr_test.imbal <- sapply(1:nreps.null, function(i) {
    sim.low <- cb.sims.sim_sigmoid(n=n, unbalancedness = 1.2, eff_sz=0)
    
    res <- cb.detect.caus_cdcorr(sim.low$Ys, sim.low$Ts, sim.low$Xs,
                                 num.threads = parallel::detectCores(),
                                 R=R)
    res$Test$p.value < alpha - eps
  })
  expect_true(mean(pwr_test.imbal) < alpha + eps)
})

test_that("throws error when no samples are retained", {
  sim.low <- cb.sims.sim_linear(n=n, eff_sz=eff_sz, unbalancedness = 10)
  
  expect_error(suppressWarnings(cb.detect.caus_cdcorr(sim.low$Ys, sim.low$Ts, sim.low$Xs,
                                                      num.threads = ncores, R=R)))
})

test_that("raises warning when few samples are retained", {
  sim.low <- cb.sims.sim_linear(n=n, eff_sz=eff_sz, unbalancedness = 3)
  
  expect_warning(cb.detect.caus_cdcorr(sim.low$Ys, sim.low$Ts, sim.low$Xs,
                                       num.threads = ncores, R=R,
                                       retain.ratio = 0.7))
})