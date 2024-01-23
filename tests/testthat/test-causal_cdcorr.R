alpha = 0.2
eps = 0.3
ncores = ifelse(testthat:::on_cran(), 1, parallel::detectCores())

test_that("reject null when balanced and ATE", {
  sim.high <- cb.sims.sim_linear(n=100, eff_sz=2)
  
  res <- cb.detect.caus_cdcorr(sim.high$Ys, sim.high$Ts, sim.high$Xs,
                               num.threads = ncores,
                               R=50)
  expect_true(res$Test$p.value < alpha)
})

test_that("reject null when imbalanced and ATE", {
  sim.low <- cb.sims.sim_linear(n=100, eff_sz=2, unbalancedness=1.5)
  
  res <- cb.detect.caus_cdcorr(sim.low$Ys, sim.low$Ts, sim.low$Xs,
                               num.threads=ncores, R=50)
  expect_true(res$Test$p.value < alpha)
})

test_that("typically fails to reject null when balanced and ATE=0", {
  skip_on_cran()
  nrep <- 50
  pwr_test <- sapply(1:nrep, function(i) {
    sim.high <- cb.sims.sim_linear(n=100, null=TRUE)
    
    res <- cb.detect.caus_cdcorr(sim.high$Ys, sim.high$Ts, sim.high$Xs,
                                 num.threads = ncores, R=50)
    res$Test$p.value < alpha
  })
  
  # should falsely reject null in favor of alternative at a rate ~alpha + small 
  # margin of error
  expect_true(mean(pwr_test) <= alpha + eps)
})

test_that("typically fails to reject null when imbalanced and ATE=0", {
  skip_on_cran()
  nrep <- 50
  pwr_test <- sapply(1:nrep, function(i) {
    sim.low <- cb.sims.sim_linear(n=100, eff_sz=2, unbalancedness = 1.5, null=TRUE)
    
    res <- cb.detect.caus_cdcorr(sim.low$Ys, sim.low$Ts, sim.low$Xs,
                                 num.threads = ncores, R=50)
    res$Test$p.value < alpha
  })
  
  expect_true(mean(pwr_test) <= alpha + eps)
})

test_that("works when covariates are given as a vector", {
  sim.high <- cb.sims.sim_linear(n=100, eff_sz=2)
  
  res <- cb.detect.caus_cdcorr(sim.high$Ys, sim.high$Ts, as.vector(sim.high$Xs), R=50)
  expect_true(res$Test$p.value < alpha)
})

test_that("works with >1 covariate levels", {
  sim.high <- cb.sims.sim_linear(n=100, eff_sz=2)
  sim.high$Xs <- cbind(sim.high$Xs, runif(100))
  
  res <- cb.detect.caus_cdcorr(sim.high$Ys, sim.high$Ts, sim.high$Xs,
                               num.threads = ncores, R=50)
  expect_true(res$Test$p.value < alpha)
})

test_that("rejects null with non-linearities and ATE", {
  sim.high <- cb.sims.sim_sigmoid(n=100, eff_sz=2)
  
  res <- cb.detect.caus_cdcorr(sim.high$Ys, sim.high$Ts, sim.high$Xs,
                               num.threads = ncores, R=50)
  expect_true(res$Test$p.value < alpha)
  
  sim.low <- cb.sims.sim_sigmoid(n=100, eff_sz=2, unbalancedness=1.5)
  
  res <- cb.detect.caus_cdcorr(sim.low$Ys, sim.low$Ts, sim.low$Xs,
                               num.threads = ncores, R=50)
  expect_true(res$Test$p.value < alpha)
})

test_that("fails to reject null with non-linearities and no ATE", {
  skip_on_cran()
  nrep <- 50
  # balanced
  pwr_test.bal <- sapply(1:nrep, function(i) {
    sim.high <- cb.sims.sim_sigmoid(n=100, null=TRUE)
    
    res <- cb.detect.caus_cdcorr(sim.high$Ys, sim.high$Ts, sim.high$Xs,
                                 num.threads = ncores, R=50)
    res$Test$p.value < alpha
  })
  expect_true(mean(pwr_test.bal) <= alpha + eps)
  
  # imbalanced
  pwr_test.imbal <- sapply(1:nrep, function(i) {
    sim.low <- cb.sims.sim_sigmoid(n=100, unbalancedness = 1.2, null=TRUE)
    
    res <- cb.detect.caus_cdcorr(sim.low$Ys, sim.low$Ts, sim.low$Xs,
                                 num.threads = parallel::detectCores(),
                                 R=50)
    res$Test$p.value < alpha
  })
  expect_true(mean(pwr_test.imbal) <= alpha + eps)
})

test_that("throws error when no samples are retained", {
  sim.low <- cb.sims.sim_linear(n=100, eff_sz=2, unbalancedness = 10)
  
  expect_error(suppressWarnings(cb.detect.caus_cdcorr(sim.low$Ys, sim.low$Ts, sim.low$Xs,
                                                      num.threads = ncores, R=50)))
})

test_that("raises warning when few samples are retained", {
  sim.low <- cb.sims.sim_linear(n=100, eff_sz=2, unbalancedness = 3)
  
  expect_warning(cb.detect.caus_cdcorr(sim.low$Ys, sim.low$Ts, sim.low$Xs,
                                       num.threads = ncores, R=50,
                                       retain.ratio = 0.7))
})