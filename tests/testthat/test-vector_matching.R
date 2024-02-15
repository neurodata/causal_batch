require(ks)

test_that("one-hot encoding works with numeric vector", {
  Ts <- c(1, 1, 0, 0)
  Ts.ohe <- as.matrix(ohe(Ts))
  expect_true(all(Ts.ohe[1,] == Ts.ohe[2,]))
  expect_true(all(Ts.ohe[3,] == Ts.ohe[4,]))
  expect_false(all(Ts.ohe[1,] == Ts.ohe[3,]))
})

test_that("one-hot encoding works with character vector", {
  Ts <- c("a", "a", "b", "b")
  Ts.ohe <- as.matrix(ohe(Ts))
  expect_true(all(Ts.ohe[1,] == Ts.ohe[2,]))
  expect_true(all(Ts.ohe[3,] == Ts.ohe[4,]))
  expect_false(all(Ts.ohe[1,] == Ts.ohe[3,]))
})

test_that("one-hot encoding works with string vector", {
  Ts <- c("aaa", "aaa", "bbb", "bbb")
  Ts.ohe <- as.matrix(ohe(Ts))
  expect_true(all(Ts.ohe[1,] == Ts.ohe[2,]))
  expect_true(all(Ts.ohe[3,] == Ts.ohe[4,]))
  expect_false(all(Ts.ohe[1,] == Ts.ohe[3,]))
})

test_that("one-hot encoding works with factor vector", {
  Ts <- factor(c("aaa", "aaa", "bbb", "bbb"))
  Ts.ohe <- as.matrix(ohe(Ts))
  expect_true(all(Ts.ohe[1,] == Ts.ohe[2,]))
  expect_true(all(Ts.ohe[3,] == Ts.ohe[4,]))
  expect_false(all(Ts.ohe[1,] == Ts.ohe[3,]))
})

test_that("zero-one distance computation works", {
  Ts <- c(1, 1, 0, 0)
  DT <- as.matrix(zero_one_dist(Ts))
  expect_true(all(DT[1:2, 1:2] == 0))
  expect_true(all(DT[3:4, 3:4] == 0))
  expect_true(all(DT[1:2, 3:4] == 1))
  expect_true(all(DT[3:4, 1:2] == 1))
})

test_that("zero-one distance computation works with arbitrary ordering and string", {
  Ts <- c("aaa", "bbb", "aaa", "bbb", "ccc")
  DT <- as.matrix(zero_one_dist(Ts))
  DT.true <- cbind(c(0, 1, 0, 1, 1), c(1, 0, 1, 0, 1), c(0, 1, 0, 1, 1),
                   c(1, 0, 1, 0, 1), c(1, 1, 1, 1, 0))
  expect_true(all(DT == DT.true))
})

test_that("vector matching with one odd-ball per-group", {
  nrep <- 100
  
  res <- sapply(1:nrep, function(i) {
    Ts <- c(rep(1, 100), rep(2, 100))
    Xs <- cbind(c(-4, runif(198), 5), c(-4, runif(198), 5), runif(200))
    
    retained.ids <- suppressWarnings(cb.align.vm_trim(Ts, Xs))
    
    excl_samps.s1 <- !(1 %in% retained.ids)
    excl_samps.s200 <- !(200 %in% retained.ids)
    
    incl_samps <- sum(2:199 %in% retained.ids)/198 > .8
    # want to exclude samples 1 and 200 and include all other samples
    # at a high rate
    return(excl_samps.s1 + excl_samps.s200 + incl_samps == 3)
  })
  # check that works most of time
  expect_true(mean(res) > .8)
})

test_that("as unbalancedness increases, fewer samples retained by VM", {
  nrep <- 20
  res <- sapply(1:nrep, function(i) {
    sim.high <- cb.sims.sim_sigmoid(n=300, unbalancedness=1)
    retained.high <- cb.align.vm_trim(sim.high$Ts, sim.high$Xs)
    
    sim.mod <- cb.sims.sim_sigmoid(n=300, unbalancedness=2)
    retained.mod <- cb.align.vm_trim(sim.mod$Ts, sim.mod$Xs)
    
    sim.low <- cb.sims.sim_sigmoid(n=300, unbalancedness=3)
    retained.low <- cb.align.vm_trim(sim.low$Ts, sim.low$Xs, retain.ratio = 0)
    
    rank.lengths <- rank(c(length(retained.high), length(retained.mod), length(retained.low)))
    return(all(rank.lengths == c(3, 2, 1)))
  })
  expect_true(mean(res) > .8)
})

test_that("VM throws warning when samples retained is low", {
  sim.low <- cb.sims.sim_sigmoid(n=200, unbalancedness=3)
  expect_warning(cb.align.vm_trim(sim.low$Ts, sim.low$Xs, retain.ratio = 0.7))
})

test_that("VM throws error when no samples retained", {
  sim.low <- cb.sims.sim_sigmoid(unbalancedness=10)
  expect_error(cb.align.vm_trim(sim.low$Ts, sim.low$Xs, retain.ratio = 0))
})

approx.overlap <- function(X1, X2, nbreaks=100) {
  xbreaks <- seq(from=-1, to=1, length.out=nbreaks)
  x1.dens <- kde(X1, eval.points=xbreaks)$estimate
  x2.dens <- kde(X2, eval.points=xbreaks)$estimate
  
  sum(pmin(x1.dens/sum(x1.dens), x2.dens/sum(x2.dens)))
}

test_that("VM increases covariate overlap", {
  sim.mod <- cb.sims.sim_sigmoid(n=200, unbalancedness = 2)
  retained.ids <- cb.align.vm_trim(sim.mod$Ts, sim.mod$Xs, retain.ratio = 0.2)
  
  Ts.tilde <- sim.mod$Ts[retained.ids]
  Xs.tilde <- sim.mod$Xs[retained.ids]
  
  ov.before <- approx.overlap(sim.mod$Xs[sim.mod$Ts == 0], sim.mod$Xs[sim.mod$Ts == 1])
  ov.after <- approx.overlap(Xs.tilde[Ts.tilde == 0], Xs.tilde[Ts.tilde == 1])
  
  expect_true(ov.before < ov.after)
})