# docker run -ti --entrypoint /bin/bash -v /cis/home/ebridge2/Documents/research/graphstats_repos/causal_batch/:/causal_batch neurodata/cdcorr:0.0.1
# cd causal_batch/docs/causal_cdcorr_paper/Figure_reproduction/
require(causalBatch)
require(ggplot2)
require(tidyverse)
require(parallel)
source("../helpers/test_statistics.R")
mode = "Deploy"
#simnames_to_run = c("Heteroskedastic", "K-Class")
simnames_to_run = c("Sigmoidal", "Non-Monotone")

sims = list("Heteroskedastic"=list(sim_fn=cb.sims.cate.heteroskedastic_sim, sim_opts=list()),
            "K-Class"=list(sim_fn=cb.sims.cate.kclass_sigmoidal_sim, sim_opts=list()), 
            "Sigmoidal"=list(sim_fn=cb.sims.cate.sigmoidal_sim, sim_opts=list()),
            "Non-Monotone"=list(sim_fn=cb.sims.cate.nonmonotone_sim, sim_opts=list()))

if(mode == "Test") {
  nrep = 5
  nbreaks = 5
  R = 100
} else if (mode == "Deploy") {
  nrep = 100
  nbreaks = 8
  R = 1000
}
n = 100
dimensionalities = list(10, 101)

batch.tests <- list("Causal cDCorr"=test.caus_cdcorr, "cDCorr"=test.cdcorr, "MANOVA"=test.manova, 
                    "GCM"=test.gcm, "RCiT"=test.rcit, "RCoT"=test.rcot, "PERMANOVA"=test.permanova,
                    "KCD"=test.kcd)

## ------- Power Simulations ------------- ##
eff_szs = seq(0, 1, length.out = nbreaks)
balances = list("High"=0.8, "Low"=0.4)
sim.tabl.power <- do.call(rbind, lapply(simnames_to_run, function(simname) {
  do.call(rbind, lapply(names(balances), function(balancename) {
    do.call(rbind, lapply(dimensionalities, function(dimen) {
      do.call(rbind, lapply(eff_szs, function(eff_sz) {
        do.call(rbind, lapply(1:nrep, function(i) {
          return(data.frame(Setting=simname, Balance=balancename, Dimensionality=dimen, Effect.Size=eff_sz, Index=i))
        }))
      }))
    }))
  }))
}))

sim.list.power <- split(sim.tabl.power, seq(nrow(sim.tabl.power)))

results.power <- do.call(rbind, lapply(sim.list.power, function(row) {
  time.st = Sys.time()
  print(row)
  sim.opts <- sims[[row$Setting]]
  sim.opts$sim_opts$n = n; sim.opts$sim_opts$p = row$Dimensionality; sim.opts$sim_opts$eff_sz = row$Effect.Size
  sim.opts$sim_opts$balance = balances[[row$Balance]]
  
  sim <- do.call(sim.opts$sim_fn, sim.opts$sim_opts)
  test.outs <- do.call(rbind, lapply(names(batch.tests), function(testname) {
    res <- tryCatch({
      do.call(batch.tests[[testname]], list(Ys=sim$Ys, Ts=sim$Ts, Xs=sim$Xs, R=R, ncores=detectCores() - 1))
    }, error=function(e) {return(list(Estimate=NaN, p.value=NaN))})
    return(data.frame(Setting=row$Setting, Balance=row$Balance, Dimensionality=row$Dimensionality,
                      Effect.Size=row$Effect.Size, Index=row$Index, Method=testname,
                      Estimate=res$Estimate, p.value=res$p.value))
  }))
  time.end = Sys.time()
  print(time.end - time.st)
  return(test.outs)
})) #, mc.cores=detectCores() - 1))

## ------- Null Simulations ----------- ##
balances = seq(1, 0.2, length.out = nbreaks)

sim.tabl.null <- do.call(rbind, lapply(names(sims), function(simname) {
  do.call(rbind, lapply(dimensionalities, function(dimen) {
    do.call(rbind, lapply(balances, function(balance) {
      do.call(rbind, lapply(1:nrep, function(i) {
        return(data.frame(Setting=simname, Dimensionality=dimen, Balance=balance, Index=i))
      }))
    }))
  }))
}))

sim.list.null <- split(sim.tabl.null, seq(nrow(sim.tabl.null)))
results.null <- do.call(rbind, mclapply(sim.list.null, function(row) {
  print(row)
  sim.opts <- sims[[row$Setting]]
  sim.opts$sim_opts$n = n; sim.opts$sim_opts$p = row$Dimensionality; sim.opts$sim_opts$balance = row$Balance
  sim.opts$sim_opts$eff_sz = 0  # specify null simulation
  
  sim <- do.call(sim.opts$sim_fn, sim.opts$sim_opts)
  test.outs <- do.call(rbind, lapply(names(batch.tests), function(testname) {
    res <- tryCatch({
      do.call(batch.tests[[testname]], list(Ys=sim$Ys, Ts=sim$Ts, Xs=sim$Xs, R=R, ncores=detectCores() - 1))
    }, error=function(e) {return(list(Estimate=NaN, p.value=NaN))})
    return(data.frame(Setting=row$Setting, Balance=row$Balance, Dimensionality=row$Dimensionality,
                      Index=row$Index, Method=testname,
                      Estimate=res$Estimate, p.value=res$p.value))
  }))
  gc()
  return(test.outs)
}, mc.cores=detectCores() - 1))