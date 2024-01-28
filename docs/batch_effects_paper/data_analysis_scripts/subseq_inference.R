# docker run -ti --entrypoint /bin/bash -v /cis/project/ndmg/batch_effects/:/data -v /cis/home/ebridge2/Documents/research/graphstats_repos/causal_batch/:/base neurodata/causal_batch:0.0.1
# docker run -ti --entrypoint /bin/bash -v /mnt/nfs2/batch_effects/:/data -v /home/eric/Documents/research/graphstats-repos/causal_batch/:/base neurodata/causal_batch:0.0.1
require(tidyverse)
require(parallel)
require(causalBatch)
require(cdcsis)
require(energy)
require(parallelDist)
require(igraph)
require(gcm)
source("../utilities/cond_alg_helpers.R")

select <- dplyr::select
in.path <- '/data/'
n.vertices <- 116
pheno.name <- "CoRR_AggregatedPhenotypicData"
pheno.path <- file.path(in.path, sprintf('phenotypic/%s.csv', pheno.name))
ncores <- detectCores() - 1
parcellation <- "AAL"
modality <- "fMRI"
cohort <- "CoRR"
am.clique <- c("NYU2", "IBATRT", "MRN1", "UWM", "NYU1")
ncores <- parallel::detectCores() - 1
#as.clique <- c("SWU4", "HNU1", "BNU3", "SWU1", "BNU2", "IPCAS1", "BNU1", "IPCAS6", "IPCAS3", "SWU2", "IPCAS4")

mri.path <- file.path(in.path, modality, parcellation)

datasets=c("UWM", "NYU1", "Utah1", "MRN1", "IBATRT", "UPSM1", "NYU2",
           "BMB1", "IPCAS8", "IPCAS4", "SWU3", "SWU2", "IACAS1", "JHNU",
           "IPCAS5", "IPCAS2", "IPCAS3", "BNU1", "IPCAS1", "BNU2",
           "SWU1", "IPCAS7", "BNU3", "XHCUMS", "HNU1", "SWU4", "NKI24tr645", "NKI24tr1400",
           "NKI24tr2500")

gr.names <- list.files(path=mri.path, pattern="*.csv", recursive=TRUE)
gr.names <- gr.names[sapply(gr.names, function(name) any(sapply(datasets, function(dataset) {str_detect(name, dataset)})))]
vertices <- 1:n.vertices
fmt <- 'ncol'

list2array <- function(x) {
  good.ids <- !sapply(x, function(xi) is.null(xi))
  x <- x[good.ids]
  return(list(good.ids=good.ids, result=t(simplify2array(x))))
}

read.gr <- function(name, mri.path='', format='ncol') {
  tryCatch({
    g <- igraph::read_graph(file.path(mri.path, name), format=format)
    V.incl <- as.numeric(V(g)$name)
    V.notincl <- (1:116)[!sapply(1:116, function(i) i %in% V.incl)]
    if (length(V.notincl) > 0) {
      g <- add_vertices(g, length(V.notincl), attr=list(name=as.character(V.notincl)))
    }
    g.perm <- permute.vertices(g, as.numeric(V(g)$name))
    g.adj <- get.adjacency(g.perm, type="both", attr="weight", sparse=FALSE)
    diag(g.adj) <- 0
    return(as.vector(g.adj))
  }, error=function(e) {
    return(NULL)
  })
}
gr.out <- list2array(mclapply(gr.names, function(name) {
  read.gr(name,  mri.path=mri.path, format=fmt)
}, mc.cores=ncores))
gr.dat.full <- gr.out$result; good.ids <- gr.out$good.ids

cov.full <- read.csv(pheno.path)
spl.names <- strsplit(basename(gr.names[good.ids]), '_|-')
dset.names <- strsplit(gr.names, '/')
cov.dat <- do.call(rbind, mclapply(1:length(spl.names), function(i) {
  spl.name <- spl.names[[i]]; dset.name <- dset.names[[i]]
  subid <- as.integer(spl.name[2]); sesid <- as.integer(spl.name[4])
  dset <- dset.name[length(dset.name)-1]
  cov.sc=cov.full %>%
    mutate(ID=row_number()) %>%
    filter(SUBID == subid) %>%
    mutate(Invalid.Entries=as.numeric(as.character(SEX) == "#") +
             as.numeric(as.character(AGE_AT_SCAN_1) == "#")) %>%
    ungroup() %>%
    filter(Invalid.Entries == min(Invalid.Entries)) %>%
    filter(ID==min(ID)) %>%
    dplyr::select(SEX, AGE_AT_SCAN_1) %>%
    rename(Sex=SEX, Age=AGE_AT_SCAN_1) %>%
    mutate(Subid=subid, Session=sesid, Dataset=dset,
           Sex=as.numeric(as.character(Sex)), Age=as.numeric(as.character(Age)))
  if (dim(cov.sc)[1] == 0) {
    return(data.frame(Subid=subid, Session=sesid, Dataset=dset, Sex=NA, Age=NA))
  } else {
    return(cov.sc)
  }
}, mc.cores=ncores))

continent <- c("IBATRT"="North America", "Utah1"="North America", "IPCAS2"="Asia", "SWU1"="Asia", "UWM"="North America", "XHCUMS"="Asia", "SWU4"="Asia",
               "BNU2"="Asia", "IPCAS3"="Asia", "SWU3"="Asia", "IPCAS4"="Asia", "NYU2"="North America", "IPCAS1"="Asia",
               "IPCAS7"="Asia", "UPSM1"="North America", "IACAS1"="Asia", "IPCAS5"="Asia", "NYU1"="North America", "NYU2"="North America", "BNU1"="Asia",
               "MRN1"="North America", "BNU3"="Asia", "HNU1"="Asia", "SWU2"="Asia", "IPCAS8"="Asia", "JHNU"="Asia", "IPCAS6"="Asia",
               "BMB1"="Europe", "NKI24tr645"="North America", "NKI24tr1400"="North America", "NKI24tr2500"="North America")
cov.dat$Continent <- as.character(continent[as.character(cov.dat$Dataset)])

cov.dat <- cov.dat %>%
  mutate(Dataset=sub("_", "", Dataset))

# strip entries with no phenotypic data
retain.ids <- complete.cases(cov.dat)
cov.dat <- cov.dat[retain.ids,] %>%
  ungroup() %>% mutate(Continent=factor(Continent),
                       Sex=factor(Sex))
gr.dat <- gr.dat.full[retain.ids,]

R=1000

fns <- list("Causal cComBat"=cb.correct.caus_cComBat, "cComBat"=cond.combat, "ComBat"=assoc.combat,
            "Raw"=raw.preproc)

covars.tbl <- cov.dat %>% ungroup() %>% mutate(id=row_number())

pos2coord<-function(pos=NULL, coord=NULL, dim.mat=NULL){
  if(is.null(pos) & is.null(coord) | is.null(dim.mat)){
    stop("must supply either 'pos' or 'coord', and 'dim.mat'")
  }
  if(is.null(pos) & !is.null(coord) & !is.null(dim.mat)){
    pos <- ((coord[,2]-1)*dim.mat[1])+coord[,1] 
    return(pos)
  }
  if(!is.null(pos) & is.null(coord) & !is.null(dim.mat)){
    coord <- matrix(NA, nrow=length(pos), ncol=2)
    coord[,1] <- ((pos-1) %% dim.mat[1]) +1
    coord[,2] <- ((pos-1) %/% dim.mat[1]) +1
    return(coord)
  }
}

cohorts=list("All"=datasets,
             "American Clique"=c("NYU2", "IBATRT", "MRN1", "UWM", "NYU1"))

# AAL parcellation has 116 vertices
nv <- 116

preproc.dat <- lapply(names(cohorts), function(crt) {
  result <- lapply(names(fns), function(fn.name) {
    tryCatch({
      subs.in.crt <- (covars.tbl %>%
                        filter(Dataset %in% cohorts[[crt]]))$id
      Ys.crt <- gr.dat[subs.in.crt,]; Ts.crt <- covars.tbl$Dataset[subs.in.crt]
      Xs.crt <- covars.tbl[subs.in.crt,] %>% select(Age, Sex, Continent) %>% mutate(Sex=factor(Sex))
      
      # if american clique, continent will be a redundant variable and cause an error
      # since it has only 1 level making the model matrix
      if (crt == "American Clique") {
        match.form <- "Age + Sex"
      } else {
        match.form <- "Age + Sex + Continent"; exact=as.formula("~Sex + Continent")
      }
      
      # "correct" batch effect
      cor.dat <- do.call(fns[[fn.name]], list(Ys.crt, Ts.crt, Xs.crt, match.form=match.form,
                                              match.args=list(method="nearest", exact=exact, 
                                                              replace=FALSE, caliper=.1)))
      
      list(Ys=cor.dat$Ys.corrected, Ts=cor.dat$Ts, Xs=cor.dat$Xs)
    }, error=function(e) {
      NULL
    })
  })
  names(result) <- names(fns)
  result
})
names(preproc.dat) <- names(cohorts)

saveRDS(preproc.dat, "../data/corrected_data.rds")

preproc.dat <- readRDS("../data/corrected_data.rds")
output <- lapply(names(cohorts), function(crt) {
  result <- lapply(names(fns), function(fn.name) {
    cor.dat <- preproc.dat[[crt]][[fn.name]]
    # test whether evidence to reject that edge and sex | age are independent
    X.sex <- cor.dat$Xs$Sex; X.age <- cor.dat$Xs$Age
    d <- 116^2
    test = do.call(rbind, mclapply(1:d, function(i) {
      i.coord <- pos2coord(i, dim.mat=c(nv, nv))
      tryCatch({
        # check if coordinate in upper triangle
        if (i.coord[1] > i.coord[2]) {
          test <- gcm(cor.dat$Ys[,i,drop=FALSE], as.matrix(causalBatch:::ohe(X.sex)), matrix(X.age, ncol=1), R=R, regr.method="xgboost")
          if(i %% 100 == 0) {
            print(i)
          }
          data.frame(Row=i.coord[1], Column=i.coord[2], Edge=i, Statistic=test$statistic,
                     p.value=test$p.value, Cohort=crt, Method=fn.name)
        } else {
          NULL
        }
      }, error=function(e) {
        data.frame(Row=i.coord[1], Column=i.coord[2], Edge=i, Statistic=NA,
                   p.value=NA, Cohort=crt, Method=fn.name)
      })
    }, mc.cores = ncores - 1))
  })
  names(result) <- names(fns)
  do.call(rbind, result)
})
output.subseq <- do.call(rbind, output)

saveRDS(outputs.subseq, "../data/subseq.rds")