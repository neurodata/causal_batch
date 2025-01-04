test_that("Matching cComBat makes outcomes for each group more similar", {
  sim <- cb.sims.sim_sigmoid(n=200, unbalancedness=1.5)
  
  xbreaks <- seq(from=-1, to=1, length.out=100)
  
  lo.before.0 <- loess(sim$Ys[sim$Ts == 0,1] ~ sim$Xs[sim$Ts == 0])
  lo.before.1 <- loess(sim$Ys[sim$Ts == 1,1] ~ sim$Xs[sim$Ts == 1])
  
  dif.before <- mean(abs(predict(lo.before.0, xbreaks) - predict(lo.before.1, xbreaks)), na.rm=TRUE)
  
  cor.sim <- cb.correct.matching_cComBat(sim$Ys, sim$Ts, data.frame(Covar=sim$Xs), 
                                     match.form="Covar")
  
  lo.after.0 <- loess(cor.sim$Ys.corrected[cor.sim$Ts == 0,1] ~ cor.sim$Xs$Covar[cor.sim$Ts == 0])
  lo.after.1 <- loess(cor.sim$Ys.corrected[cor.sim$Ts == 1,1] ~ cor.sim$Xs$Covar[cor.sim$Ts == 1])
  
  dif.after <- mean(abs(predict(lo.after.0, xbreaks) - predict(lo.after.1, xbreaks)), na.rm=TRUE)
  
  expect_true(dif.before > dif.after)
})

test_that("Matching cComBat raises warning when few samples are retained", {
  sim.low <- cb.sims.sim_linear(n=200, unbalancedness = 2.5)
  
  expect_warning(cb.correct.matching_cComBat(sim.low$Ys, sim.low$Ts, data.frame(Covar=sim.low$Xs), 
                                         match.form="Covar", retain.ratio=0.8))
})

test_that("Matching cComBat raises error when no samples are retained", {
  sim.low <- cb.sims.sim_linear(n=100, unbalancedness = 10)
  
  expect_error(suppressWarnings(cb.correct.matching_cComBat(sim.low$Ys, sim.low$Ts, data.frame(Covar=sim.low$Xs), 
                                         match.form="Covar", retain.ratio=0.5)))
})

test_that("apply Matching cComBat raises error when oos data has inconsistent batches", {
  sim <- cb.sims.sim_sigmoid(n=200, unbalancedness=1.5)
  
  fit_ccombat <- cb.correct.matching_cComBat(sim$Ys, sim$Ts, data.frame(Covar=sim$Xs),
                                         match.form="Covar")
  
  Ys <- rbind(sim$Ys, c(1, 0)); Ts <- c(sim$Ts, 2); Xs <- c(sim$Xs, 1)
  expect_error(cb.correct.apply_cComBat(Ys, Ts, data.frame(Covar=sim$Xs), Model=fit_ccombat$Model))
})

test_that("AIPW cComBat makes outcomes for each group more similar", {
  sim <- cb.sims.sim_sigmoid(n=200, unbalancedness=1.5)
  
  xbreaks <- seq(from=-1, to=1, length.out=100)
  
  lo.before.0 <- loess(sim$Ys[sim$Ts == 0,1] ~ sim$Xs[sim$Ts == 0])
  lo.before.1 <- loess(sim$Ys[sim$Ts == 1,1] ~ sim$Xs[sim$Ts == 1])
  
  dif.before <- mean(abs(predict(lo.before.0, xbreaks) - predict(lo.before.1, xbreaks)), na.rm=TRUE)
  
  cor.sim <- cb.correct.aipw_cComBat(sim$Ys, sim$Ts, data.frame(Covar=sim$Xs), 
                                     "Covar")
  
  lo.after.0 <- loess(cor.sim$Ys.corrected[cor.sim$Ts == 0,1] ~ cor.sim$Xs$Covar[cor.sim$Ts == 0])
  lo.after.1 <- loess(cor.sim$Ys.corrected[cor.sim$Ts == 1,1] ~ cor.sim$Xs$Covar[cor.sim$Ts == 1])
  
  dif.after <- mean(abs(predict(lo.after.0, xbreaks) - predict(lo.after.1, xbreaks)), na.rm=TRUE)
  
  expect_true(dif.before > dif.after)
})

test_that("AIPW cComBat raises warning when few samples are retained", {
  sim.low <- cb.sims.sim_linear(n=200, unbalancedness = 2.5)
  
  expect_warning(cb.correct.aipw_cComBat(sim.low$Ys, sim.low$Ts, data.frame(Covar=sim.low$Xs), 
                                         aipw.form="Covar", retain.ratio=0.8))
})

test_that("AIPW cComBat raises error when no samples are retained", {
  sim.low <- cb.sims.sim_linear(n=100, unbalancedness = 10)
  
  expect_error(suppressWarnings(cb.correct.aipw_cComBat(sim.low$Ys, sim.low$Ts, data.frame(Covar=sim.low$Xs), 
                                                        aipw.form="Covar", retain.ratio=0.5)))
})

test_that("apply AIPW cComBat works reasonably well", {
  sim <- cb.sims.sim_sigmoid(n=200, unbalancedness=1.5)
  
  xbreaks <- seq(from=-1, to=1, length.out=100)
  
  lo.before.0 <- loess(sim$Ys[sim$Ts == 0,1] ~ sim$Xs[sim$Ts == 0])
  lo.before.1 <- loess(sim$Ys[sim$Ts == 1,1] ~ sim$Xs[sim$Ts == 1])
  
  dif.before <- mean(abs(predict(lo.before.0, xbreaks) - predict(lo.before.1, xbreaks)), na.rm=TRUE)
  
  cor.sim <- cb.correct.aipw_cComBat(sim$Ys, sim$Ts, data.frame(Covar=sim$Xs), 
                                     "Covar")
  
  applied.cor.sim <- cb.correct.apply_aipw_cComBat(sim$Ys, sim$Ts, data.frame(Covar=sim$Xs), cor.sim$Model)
  
  lo.after.0 <- loess(applied.cor.sim[sim$Ts == 0,1] ~ sim$Xs[sim$Ts == 0])
  lo.after.1 <- loess(applied.cor.sim[sim$Ts == 1,1] ~ sim$Xs[sim$Ts == 1])
  
  dif.after <- mean(abs(predict(lo.after.0, xbreaks) - predict(lo.after.1, xbreaks)), na.rm=TRUE)
  
  expect_true(dif.before > dif.after)
})

test_that("apply AIPW cComBat raises error when new data has inconsistent batches", {
  sim <- cb.sims.sim_sigmoid(n=200, unbalancedness=1.5)
  
  fit_ccombat <- cb.correct.aipw_cComBat(sim$Ys, sim$Ts, data.frame(Covar=sim$Xs),
                                             aipw.form="Covar")
  
  Ys <- rbind(sim$Ys, c(1, 0)); Ts <- c(sim$Ts, 2); Xs <- c(sim$Xs, 1)
  expect_error(cb.correct.apply_aipw_cComBat(Ys, Ts, data.frame(Covar=sim$Xs), fit_ccombat$Model))
})