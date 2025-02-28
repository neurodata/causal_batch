# docker run -ti --entrypoint /bin/bash -v /cis/home/ebridge2/Documents/research/graphstats_repos/causal_batch/:/causal_batch neurodata/cdcorr:0.0.1
# cd causal_batch/docs/causal_cdcorr_paper/Figure_reproduction/
source("../helpers/test_statistics.R")
require(causalBatch)
require(tidyverse)
require(energy)
require(jsonlite)
require(np)
require(parallel)
modality = "dMRI"
null.splits <- 5  # number of null splits for each dataset

dcorr_ksample <- function(Ts, Xs) {
  # ts is the distance matrix for the grouping variable (here, batches)
  # xs is the covariate matrix, which should already be one-hot encoded
  
  # compute zero-one distance for the categorical group assignment
  Ts <- as.factor(Ts)
  DT <- as.matrix(causalBatch:::zero_one_dist(Ts, levels=levels(Ts))$DT)
  
  DX <- as.matrix(dist(as.matrix(Xs), method = "euclidean"))
  
  dcor(DX, DT)
}

inputs.modality <- readRDS("../data/abcd/preproc_abcd_data.csv")[[modality]]
all.dat <- inputs.modality$All
Xs.full <- inputs.modality$Xs; Ys.full <- inputs.modality$Ys; Xs.scaled.full <- inputs.modality$Xs.scaled
Ts.full <- inputs.modality$Ts; sites <- unique(Ts.full)
bws <- inputs.modality$bws.scott

site_details <- all.dat %>%
  select(site, scanner_manufacturer, scanner_model) %>%
  group_by(site, scanner_manufacturer, scanner_model) %>%
  summarise(n_subjects = n()) %>%  # or n_distinct(participant_id) if you want to ensure unique subjects
  ungroup() %>%
  group_by(site) %>%
  arrange(site, desc(n_subjects), .by_group=TRUE)

maximal_scanner <- site_details %>%
  group_by(site) %>%
  slice_max(n_subjects, n=1)

batch.tests <- list("Causal cDCorr"=test.caus_cdcorr, "cDCorr"=test.cdcorr, "MANOVA"=test.manova, 
                    "GCM"=test.gcm, "RCIT"=test.rcit, "RCoT"=test.rcot, "PERMANOVA"=test.permanova,
                    "KCD"=test.kcd)

##------------ NULL Experiments ---------------##
# Create all site-split combinations first
null.site_split_combos <- expand.grid(
  site = sites,
  split_idx = 1:null.splits,
  stringsAsFactors = FALSE
)

results.null <- do.call(rbind, lapply(1:nrow(null.site_split_combos), function(i) {
  sitevar <- null.site_split_combos$site[i]
  split_idx <- null.site_split_combos$split_idx[i]
  
  print(sprintf("Site: %s, Split: %d", sitevar, split_idx))
  
  site_scanner_to_use <- maximal_scanner %>% dplyr::filter(site == sitevar)
  idx.site <- which(Ts.full == sitevar, 
                    all.dat$scanner_manufacturer == site_scanner_to_use$scanner_manufacturer,
                    all.dat$scanner_model == site_scanner_to_use$scanner_model)
  
  Ys <- Ys.full[idx.site,]
  Xs <- Xs.full[idx.site,]
  Xs.scaled <- Xs.scaled.full[idx.site,]
  
  # randomly split the dataset into two batches
  Ts <- sample(c(1, 2), prob=c(0.5, 0.5), replace=TRUE, size=length(idx.site))
  
  # compute covariate overlap of the split dataset
  overlap <- dcorr_ksample(Ts, Xs.scaled)
  
  do.call(rbind, lapply(names(batch.tests), function(test) {
    print(test)
    res <- tryCatch({
      do.call(batch.tests[[test]], list(Ys=Ys, Ts=Ts, Xs=Xs, normalize=TRUE, width=bws, ncores=detectCores() - 1))
    }, error=function(e) {return(list(Estimate=NaN, p.value=NaN))})
    
    return(data.frame(
      Dataset=sitevar, 
      index=split_idx, 
      Method=test, 
      Estimate=res$Estimate, 
      p.value=res$p.value, 
      Overlap.SS=overlap
    ))
  }))
}))


##------------ ALTERNATIVE Experiments ---------------##
# Create all site-split combinations first
site_pairs <- expand.grid(site1 = maximal_scanner$site, 
                          site2 = maximal_scanner$site) %>%
  # Remove self-pairs
  filter(site1 != site2) %>%
  # Remove duplicates (e.g., CHLA-CUB vs CUB-CHLA)
  filter(as.character(site1) < as.character(site2)) %>%
  left_join(maximal_scanner %>%
              rename(scanner_model_1=scanner_model, scanner_manufacturer_1=scanner_manufacturer), 
            by = c("site1" = "site")) %>%
  left_join(maximal_scanner %>%
              rename(scanner_model_2=scanner_model, scanner_manufacturer_2=scanner_manufacturer), 
            by = c("site2" = "site")) %>%
  filter(scanner_model_1 != scanner_model_2 | scanner_manufacturer_1 != scanner_manufacturer_2)

# Next create all combinations of site pairs and iterations
alt.site_split_combos <- do.call(rbind, lapply(1:nrow(site_pairs), function(i) {
  data.frame(
    pair_idx = i,
    split_idx = 1:null.splits,
    site1 = site_pairs$site1[i],
    site2 = site_pairs$site2[i],
    scanner_manufacturer_1 = site_pairs$scanner_manufacturer_1[i],
    scanner_model_1 = site_pairs$scanner_model_1[i],
    scanner_manufacturer_2 = site_pairs$scanner_manufacturer_2[i],
    scanner_model_2 = site_pairs$scanner_model_2[i]
  )
}))


results.alt <- do.call(rbind, mclapply(1:nrow(alt.site_split_combos), function(i) {
  combo <- alt.site_split_combos[i,]
  print(sprintf("Sites: %s vs %s, Split: %d", combo$site1, combo$site2, combo$split_idx))
  
  # Get indices for first site
  idx.site1 <- which(Ts.full == combo$site1 & 
                       all.dat$scanner_manufacturer == combo$scanner_manufacturer_1 &
                       all.dat$scanner_model == combo$scanner_model_1)
  
  # Get indices for second site
  idx.site2 <- which(Ts.full == combo$site2 & 
                       all.dat$scanner_manufacturer == combo$scanner_manufacturer_2 &
                       all.dat$scanner_model == combo$scanner_model_2)
  
  # Randomly select half of each site's data
  idx.site1.subset <- sample(idx.site1, size=floor(length(idx.site1)/2))
  idx.site2.subset <- sample(idx.site2, size=floor(length(idx.site2)/2))
  
  # Combine the subsetted data
  idx.combined <- c(idx.site1.subset, idx.site2.subset)
  Ts <- c(rep(1, length(idx.site1.subset)), rep(2, length(idx.site2.subset)))
  Ys <- Ys.full[idx.combined,]
  Xs <- Xs.full[idx.combined,]
  Xs.scaled <- Xs.scaled.full[idx.combined,]
  
  # compute covariate overlap
  overlap <- dcorr_ksample(Ts, Xs.scaled)
  
  do.call(rbind, lapply(names(batch.tests), function(test) {
    print(test)
    res <- tryCatch({
      do.call(batch.tests[[test]], list(Ys=Ys, Ts=Ts, Xs=Xs, width=bws, normalize=TRUE))
    }, error=function(e) {return(list(Estimate=NaN, p.value=NaN))})
    
    gc()
    return(data.frame(
      Dataset1=combo$site1,
      Dataset2=combo$site2,
      index=combo$split_idx,
      Method=test,
      Estimate=res$Estimate,
      p.value=res$p.value,
      Overlap.SS=overlap
    ))
  }))
}, mc.cores=detectCores() - 1))

saveRDS(list(Null=results.null, Alt=results.alt), "../data/abcd/results.rds")