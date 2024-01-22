test_that("Causal cComBat makes outcomes for each group more similar", {
  sim <- cb.sims.sim_sigmoid(n=200, unbalancedness=1.5)
  
  xbreaks <- seq(from=-1, to=1, length.out=100)
  
  lo.before.0 <- loess(sim$Ys[sim$Ts == 0,1] ~ sim$Xs[sim$Ts == 0])
  lo.before.1 <- loess(sim$Ys[sim$Ts == 1,1] ~ sim$Xs[sim$Ts == 1])
  
  dif.before <- mean(abs(predict(lo.before.0, xbreaks) - predict(lo.before.1, xbreaks)), na.rm=TRUE)
  
  cor.sim <- cb.correct.caus_cComBat(sim$Ys, sim$Ts, data.frame(Covar=sim$Xs), 
                                     match.form="Covar")
  
  lo.after.0 <- loess(cor.sim$Ys.corrected[cor.sim$Ts == 0,1] ~ cor.sim$Xs$Covar[cor.sim$Ts == 0])
  lo.after.1 <- loess(cor.sim$Ys.corrected[cor.sim$Ts == 1,1] ~ cor.sim$Xs$Covar[cor.sim$Ts == 1])
  
  dif.after <- mean(abs(predict(lo.after.0, xbreaks) - predict(lo.after.1, xbreaks)), na.rm=TRUE)
  
  expect_true(dif.before > dif.after)
})

test_that("raises warning when few samples are retained", {
  sim.low <- cb.sims.sim_linear(n=200, unbalancedness = 2)
  
  expect_warning(cb.correct.caus_cComBat(sim.low$Ys, sim.low$Ts, data.frame(Covar=sim.low$Xs), 
                                         match.form="Covar", retain.ratio=0.7))
})

test_that("raises error when no samples are retained", {
  sim.low <- cb.sims.sim_linear(n=100, unbalancedness = 10)
  
  expect_error(suppressWarnings(cb.correct.caus_cComBat(sim.low$Ys, sim.low$Ts, data.frame(Covar=sim.low$Xs), 
                                         match.form="Covar", retain.ratio=0.5)))
})