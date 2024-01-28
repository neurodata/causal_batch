require(matrixStats)
require(BiocParallel)
require(tidyverse)

gen.simulations <- function(sim_fn, paramlist=overlap.settings) {
  simulations=lapply(paramlist, function(params) do.call(sim_fn, params))
  names(simulations) <- names(paramlist)
  return(simulations)
}

get.scatterplot_df <- function(raw_dat) {
  do.call(rbind, lapply(names(raw_dat), function(sim_setting) {
    raw_sim <- raw_dat[[sim_setting]]
    data.frame(X=raw_sim$Xs, Y=raw_sim$Ys[,1], Batch=raw_sim$Ts, Setting=sim_setting, Overlap=raw_sim$Overlap)
  }))
}

plt.raw_dat <- function(raw_df, title="") {
  raw_df %>% 
    ggplot(aes(x=X, y=Y, color=factor(Batch))) +
    geom_point(alpha=0.5) +
    facet_grid(Setting~., switch="y") +
    scale_x_continuous(name="Covariate", breaks=c(-1, 0, 1)) +
    scale_y_continuous(name="Outcome", limits=c(min(raw_df$Y) - .1, max(raw_df$Y) + .1)) +
    scale_color_manual(name="Batch", values=batch.cols) +
    theme_bw() +
    ggtitle(title) +
    theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
          strip.background=element_blank(), text=element_text(size=13),
          axis.text.y=element_blank(), axis.ticks.y=element_blank())
}

get.true_signal <- function(raw_dat) {
  result <- lapply(names(raw_dat), function(sim_setting) {
    raw_sim <- raw_dat[[sim_setting]]
    list(Xs=raw_sim$Xtrue, Ys=raw_sim$Ytrue, Ts=raw_sim$Ttrue, Setting=sim_setting, Overlap=raw_sim$Overlap)
  })
  names(result) <- names(raw_dat)
  return(result)
}

get.true_signal_df <- function(sig_list) {
  do.call(rbind, lapply(names(sig_list), function(sim_setting) {
    sig_sim <- sig_list[[sim_setting]]
    return(data.frame(X=sig_sim$Xs, Y=sig_sim$Ys[,1], Batch=sig_sim$Ts, Setting=sig_sim$Setting, Overlap=sig_sim$Overlap))
  })) %>% mutate(Signal="Expected Signal")
}

fit.cComBat <- function(raw_dat) {
  fit_mods <- lapply(raw_dat, function(sim) {
    cor.sim.cb <- ComBat_return_model(t(sim$Ys), as.vector(sim$Ts), mod=model.matrix(~1 + sim$Xs))
    cor.sim.cb$xbounds <- c(min(sim$Xs), max(sim$Xs))
    return(cor.sim.cb)
  })
  names(fit_mods) <- names(raw_dat)
  return(fit_mods)
}

fit.ComBatGAM <- function(raw_dat) {
  fit_mods <- lapply(raw_dat, function(sim) {
    cor.sim.cb <- ComBat.GAM(sim$Ys, as.vector(sim$Ts), data.frame(Covariate=sim$Xs), 
                             nh.args = list(smooth_terms="Covariate"))
    cor.sim.cb$Model$xbounds <- c(min(sim$Xs), max(sim$Xs))
    return(cor.sim.cb)
  })
  names(fit_mods) <- names(raw_dat)
  return(fit_mods)
}

fit.caus_cComBat <- function(raw_dat) {
  fit_mods <- lapply(names(raw_dat), function(simn) {
    sim <- raw_dat[[simn]]
    tryCatch({
      retain.ids <- unique(do.call(cb.align.kway_match, list(sim$Ts, data.frame(Covar=sim$Xs), "Covar",
                                                       match.args=list(method="nearest", exact=NULL,
                                                                       replace=FALSE, caliper=.1))))
      X.tilde <- sim$Xs[retain.ids]; Y.tilde <- sim$Ys[retain.ids,]; T.tilde <- sim$Ts[retain.ids]
      
      cor.sim.ccb <- ComBat_return_model(t(Y.tilde), as.vector(T.tilde), mod=model.matrix(~1 + X.tilde))
      cor.sim.ccb$xbounds <- c(min(X.tilde), max(X.tilde))
      return(cor.sim.ccb)
    }, error=function(e) {message(sprintf("Caus. cComBat failed %s.", simn)); return(NULL)})
  })
  names(fit_mods) <- names(raw_dat)
  return(fit_mods)
}

fit.linfit_lines <- function(true_sig, mod_fits) {
  if (any(names(true_sig) != names(mod_fits))) {
    stop("You have not fit a model for the same set of simulations.")
  }
  do.call(rbind, lapply(names(true_sig), function(sim_setting) {
    sim.expected.setting <- true_sig[[sim_setting]]
    mod.setting <- mod_fits[[sim_setting]]
    retain.ids <- which(sim.expected.setting$X >= mod.setting$xbounds[1] & 
                          sim.expected.setting$X <= mod.setting$xbounds[2])
    X.tilde <- sim.expected.setting$X[retain.ids]; Y.tilde <- sim.expected.setting$Y[retain.ids,0]
    T.tilde <- sim.expected.setting$Batch[retain.ids]
    T.design <- causalBatch:::ohe(as.vector(T.tilde)); colnames(T.design) <- c(sapply(c(0, 1), function(b) {sprintf("batch%d", b)}))
    
    fit.true_covars <- (as.matrix(cbind(T.design, X.tilde)) %*% mod.setting$B.hat)[,1]
    return(data.frame(X=X.tilde, Batch=T.tilde, Y=fit.true_covars, Setting=sim_setting, Overlap=sim.expected.setting$Overlap))
  })) %>% mutate(Signal="cComBat fit")
}

plt.exp_sigs <- function(true_df, linfit_df, title="") {
  batch_effect_df <- true_df %>%
    filter(Signal == "Expected Signal") %>%
    pivot_wider(names_from=Batch, values_from=Y) %>%
    group_by(X, Setting, Overlap) %>%
    mutate(ymin=min(`0`, `1`), ymax=max(`0`, `1`))
  
  ggplot(rbind(true_df, linfit_df), aes(x=X)) +
    geom_ribbon(data=batch_effect_df, aes(ymin=ymin, ymax=ymax), alpha=0.3, color="gray") +
    geom_line(linewidth=1.2, aes(y=Y, color=factor(Batch), linetype=factor(Signal))) +
    facet_grid(Setting~., switch="y") +
    scale_linetype_manual(name="Batch", values=c(2, 1)) +
    scale_x_continuous(name="Covariate", breaks=c(-1, 0, 1)) +
    scale_y_continuous(name="Outcome", limits=c(min(sig.raw_df$Y) - .1, max(sig.raw_df$Y) + .1)) +
    scale_color_manual(name="Batch", values=batch.cols) +
    theme_bw() +
    ggtitle(title) +
    theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
          strip.background=element_blank(), text=element_text(size=13),
          axis.text.y=element_blank(), axis.ticks.y=element_blank())
}

apply.ComBat.GAM <- function(true_sig, mod_fits) {
  if (any(names(true_sig) != names(mod_fits))) {
    stop("You have not fit a model for the same set of simulations.")
  }
  do.call(rbind, lapply(names(true_sig), function(sim_setting) {
    message(sim_setting)
    sim.expected.setting <- true_sig[[sim_setting]]
    mod.setting <- mod_fits[[sim_setting]]$Model
    if (is.null(mod.setting)) {
      message(sprintf("No model fit for %s. Skipping...", sim_setting))
      return(data.frame(X=NA, Batch=NA, Y=NA, Setting=sim_setting, Overlap=sim.expected.setting$Overlap))
    }
    retain.ids <- which(sim.expected.setting$Xs >= mod.setting$xbounds[1] & 
                          sim.expected.setting$Xs <= mod.setting$xbounds[2])
    X.tilde <- sim.expected.setting$Xs[retain.ids]; Y.tilde <- sim.expected.setting$Ys[retain.ids,]
    T.tilde <- sim.expected.setting$Ts[retain.ids]
    
    fit.exp <- ComBat.GAM.apply(Y.tilde, batches=T.tilde, data.frame(Covariate=X.tilde), mod.setting)
    
    return(data.frame(X=X.tilde, Batch=T.tilde, Y=fit.exp$Data[,1], Setting=sim_setting, Overlap=sim.expected.setting$Overlap))
  })) %>% mutate(Signal="Expected Signal")
}

apply.cComBat <- function(true_sig, mod_fits) {
  if (any(names(true_sig) != names(mod_fits))) {
    stop("You have not fit a model for the same set of simulations.")
  }
  do.call(rbind, lapply(names(true_sig), function(sim_setting) {
    message(sim_setting)
    sim.expected.setting <- true_sig[[sim_setting]]
    mod.setting <- mod_fits[[sim_setting]]
    if (is.null(mod.setting)) {
      message(sprintf("No model fit for %s. Skipping...", sim_setting))
      return(data.frame(X=NA, Batch=NA, Y=NA, Setting=sim_setting, Overlap=sim.expected.setting$Overlap))
    }
    retain.ids <- which(sim.expected.setting$Xs >= mod.setting$xbounds[1] & 
                          sim.expected.setting$Xs <= mod.setting$xbounds[2])
    X.tilde <- sim.expected.setting$Xs[retain.ids]; Y.tilde <- sim.expected.setting$Ys[retain.ids,]
    T.tilde <- sim.expected.setting$Ts[retain.ids]
    T.design <- causalBatch:::ohe(as.vector(T.tilde)); colnames(T.design) <- c(sapply(c(0, 1), function(b) {sprintf("batch%d", b)}))
    
    T.tilde.locs <- lapply(c(0, 1), function(i) which(T.tilde == i))
    fit.exp <- apply_batch_adjustment(data=t(Y.tilde),
                                      batches = T.tilde.locs, design = T.design, gamma = mod.setting$Gamma, 
                                      delta = mod.setting$Delta, N=mod.setting$N, var.pooled = mod.setting$Var, 
                                      grand.mean = mod.setting$Grand.mean, B.hat = mod.setting$B.hat, 
                                      mod = model.matrix(~1 + X.tilde))
    return(data.frame(X=X.tilde, Batch=T.tilde, Y=fit.exp[,1], Setting=sim_setting, Overlap=sim.expected.setting$Overlap))
  })) %>% mutate(Signal="Expected Signal")
}

plt.corrected_true <- function(cor_true_adj, title="") {
  resid_df <- cor_true_adj %>%
    pivot_wider(names_from=Batch, values_from=Y) %>%
    group_by(X, Setting, Overlap) %>%
    mutate(ymin=min(`0`, `1`), ymax=max(`0`, `1`))
  
  ggplot(cor_true_adj, aes(x=X)) +
    geom_ribbon(data=resid_df, aes(ymin=ymin, ymax=ymax), alpha=0.3, fill="red") +
    geom_line(linewidth=1.2, aes(y=Y, color=factor(Batch))) +
    facet_grid(Setting~., switch="y") +
    scale_linetype_manual(name="Batch", values=c(2, 1)) +
    scale_x_continuous(name="Covariate", breaks=c(-1, 0, 1)) +
    scale_y_continuous(name="Outcome", limits=c(min(sig.raw_df$Y) - .1, max(sig.raw_df$Y) + .1)) +
    scale_color_manual(name="Batch", values=batch.cols) +
    theme_bw() +
    ggtitle(title) +
    theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
          strip.background=element_blank(), text=element_text(size=13),
          axis.text.y=element_blank(), axis.ticks.y=element_blank())
}

compute.signal_disparity.cate <- function(simfn, paramlist=sim.ov.settings, fname=NULL) {
  tryCatch({
    raw_dat <- gen.simulations(simfn, paramlist=paramlist)
    if (!is.null(fname)) {
      saveRDS(raw_dat, fname)
    }
    true_list <- get.true_signal(raw_dat)
    
    fit_combat <- fit.cComBat(raw_dat)
    fit_caus_combat <- fit.caus_cComBat(raw_dat)
    fit_gam_combat <- fit.ComBatGAM(raw_dat)
    
    cc_true_adj <- apply.cComBat(true_list, fit_combat)
    caus.cc_true_adj <- apply.cComBat(true_list, fit_caus_combat)
    gam.cc_true_adj <- apply.ComBat.GAM(true_list, fit_gam_combat)
    
    get.resid_avg <- function(df) {
      tryCatch({
        df %>%
          pivot_wider(names_from=Batch, values_from=Y) %>%
          group_by(X, Setting, Overlap) %>%
          mutate(Batch.Effect=`1` - `0`) %>%
          group_by(Setting, Overlap) %>%
          summarize(Avg.Effect=mean(Batch.Effect, na.rm=TRUE))
      }, error=function(e) {
        return(data.frame(Setting=cc_true_adj$Setting, Overlap=cc_true_adj$Overlap,
                          Avg.Effect=NA))
      })
    }
    cc_resid_df <- get.resid_avg(cc_true_adj) %>% mutate(Strategy="cComBat")
    caus.cc_resid_df <- get.resid_avg(caus.cc_true_adj) %>% mutate(Strategy="Causal cComBat")
    gam.cc_resid_df <- get.resid_avg(gam.cc_true_adj) %>% mutate(Strategy="cComBat-GAM")
    cc_resid_df %>% 
      rename(`cComBat`=Avg.Effect) %>%
      dplyr::select(-Strategy) %>%
      left_join(caus.cc_resid_df %>% rename(`Causal cComBat`=Avg.Effect) %>% dplyr::select(-Strategy),
                by="Setting") %>%
      left_join(gam.cc_resid_df %>% rename(`cComBat-GAM`=Avg.Effect) %>% dplyr::select(-Strategy),
                by="Setting")
  }, error=function(e) {print(e); return(list(error=raw_dat))})
}

se <- function(x) {
  x <- x[!is.nan(x)]
  sd(x)/length(x)
}

compute.mean_CATE <- function(simfn, R=1000, ncores=8, paramlist=sim.ov.settings, fname=NULL) {
  do.call(rbind, mclapply(1:R, function(r) {
    if (is.null(fname)) {
      file.name = NULL
    } else {
      file.name = sprintf("%s_%d.rds", fname, r)
    }
    
    result <- suppressWarnings(
      suppressMessages(compute.signal_disparity.cate(simfn, paramlist=paramlist, fname=file.name)
      )
    )
    if(is.data.frame(result)) {
      result <- result %>% mutate(idx=r)
    } else {
      result <- NULL
    }
    return(result)
  }, mc.cores = ncores))
}

plt.combine_plots <- function(col1, col2, col3, col4,
                              title="") {
  sep_plts <- lapply(list(col1, col2, col3, col4), function(col) {
    leg = get_legend(col)
    plt <- col + theme(legend.position="none")
    return(list(Plot=plt, Legend=leg))
  })
  
  plots <- do.call(ggarrange, c(lapply(sep_plts, function(i) {i$Plot}), ncol=length(sep_plts)))
  
  return(annotate_figure(plots, top=text_grob(title, size=20, hjust=0, x=0)))
}