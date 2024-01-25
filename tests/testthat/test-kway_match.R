require(ks)

test_that("k-way matching with one odd-ball per-group", {
  nrep <- 100
  
  res <- sapply(1:nrep, function(i) {
    Ts <- c(rep(1, 100), rep(2, 100))
    Xs <- cbind(c(-4, runif(99), runif(99), 5), runif(200), runif(200))
    
    retained.ids <- suppressWarnings(cb.align.kway_match(Ts, data.frame(Covar=Xs), match.form = "Covar.1 + Covar.2 + Covar.3",
                                                         match.args=list(method="nearest", caliper=0.5, exact=NULL,replace=FALSE)))
    
    excl_samps.s1 <- !(1 %in% retained.ids)
    excl_samps.s200 <- !(200 %in% retained.ids)
    
    incl_samps <- sum(2:199 %in% retained.ids)/198 > .9
    # want to exclude samples 1 and 200 and include all other samples
    # at a high rate
    return(excl_samps.s1 + excl_samps.s200 + incl_samps == 3)
  })
  # check that works most of time
  expect_true(mean(res) > .8)
})

test_that("as unbalancedness increases, fewer samples retained by k-way matching", {
  sim.high <- cb.sims.sim_sigmoid(n=200, unbalancedness=1)
  retained.high <- cb.align.kway_match(sim.high$Ts, data.frame(Covar=sim.high$Xs),
                                       match.form="Covar")
  
  sim.mod <- cb.sims.sim_sigmoid(n=200, unbalancedness=1.5)
  retained.mod <- cb.align.kway_match(sim.mod$Ts, data.frame(Covar=sim.mod$Xs),
                                      match.form="Covar")
  
  sim.low <- cb.sims.sim_sigmoid(n=200, unbalancedness=2.5)
  retained.low <- cb.align.kway_match(sim.low$Ts, data.frame(Covar=sim.low$Xs),
                                      match.form="Covar", retain.ratio=0)
  
  rank.lengths <- rank(c(length(retained.high), length(retained.mod), length(retained.low)))
  expect_true(all(rank.lengths == c(3, 2, 1)))
})


test_that("K-way matching throws warning when samples retained is low", {
  sim.low <- cb.sims.sim_sigmoid(n=200, unbalancedness=2)
  expect_warning(cb.align.kway_match(sim.low$Ts, data.frame(Covar=sim.low$Xs),
                                     match.form="Covar", retain.ratio = 0.7))
})

test_that("K-way matching throws error when no samples retained", {
  sim.low <- cb.sims.sim_sigmoid(unbalancedness=10)
  expect_error(suppressWarnings(cb.align.kway_match(sim.low$Ts, data.frame(Covar=sim.low$Xs),
                                                    match.form="Covar", retain.ratio = 0.5)))
})

approx.overlap <- function(X1, X2, nbreaks=100) {
  xbreaks <- seq(from=-1, to=1, length.out=nbreaks)
  x1.dens <- kde(X1, eval.points=xbreaks)$estimate
  x2.dens <- kde(X2, eval.points=xbreaks)$estimate
  
  sum(pmin(x1.dens/sum(x1.dens), x2.dens/sum(x2.dens)))
}

test_that("K-way matching increases covariate overlap", {
  sim.mod <- cb.sims.sim_sigmoid(n=200, unbalancedness = 2)
  retained.ids <- cb.align.kway_match(sim.mod$Ts, data.frame(Covar=sim.mod$Xs),
                                      match.form="Covar", retain.ratio = 0)
  
  Ts.tilde <- sim.mod$Ts[retained.ids]
  Xs.tilde <- sim.mod$Xs[retained.ids]
  
  ov.before <- approx.overlap(sim.mod$Xs[sim.mod$Ts == 0], sim.mod$Xs[sim.mod$Ts == 1])
  ov.after <- approx.overlap(Xs.tilde[Ts.tilde == 0], Xs.tilde[Ts.tilde == 1])
  
  expect_true(ov.before < ov.after)
})