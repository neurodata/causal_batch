require(tidyverse)
require(parallel)
require(causalBatch)
require(cdcsis)
require(energy)
source("./cdcorr_combat_helpers.R")

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
gr.dat.full <- gr.dat.full[retain.ids,]

retain.dims <- sapply(1:dim(gr.dat.full)[2], function(j) {
  all(sapply(unique(cov.dat$Dataset), function(dataset) {
    var(gr.dat.full[cov.dat$Dataset == dataset,j]) > 0
  }))
})
gr.dat <- gr.dat.full[,retain.dims]

R=10000

fns <- list("Adjusted (Causal)"=cb.detect.caus_cdcorr, "Conditional (Non-Causal)"=cond.dcorr, "Associational"=dcorr)

covars.tbl <- cov.dat %>% ungroup() %>% mutate(id=row_number())

Dmtx.norm <- as.matrix(parDist(gr.dat, threads=ncores))

dset.pairs <- combn(datasets)
output <- do.call(rbind, mclapply(1:dim(dset.pairs)[2], function(i) {
  tryCatch({
    dset.1 <- dset.pairs[1,x]
    dset.2 <- dset.pairs[2,x]
    print(sprintf("%d: %s, %s", x, dset.1, dset.2))
    n.1 <- sum(cov.dat$Dataset == dset.1)
    n.2 <- sum(cov.dat$Dataset == dset.2)
    if (n.1 < n.2) {
      dset.i <- dset.1; dset.j <- dset.2
      n.i <- n.1; n.j <- n.2
    } else if (n.1 >= n.2) {
      dset.i <- dset.2; dset.j <- dset.1
      n.i <- n.2; n.j <- n.1
    }
    
    cov.dat.ij <- cov.dat %>%
      filter(Dataset %in% c(dset.i, dset.j)) %>%
      mutate(Treatment = as.numeric(Dataset == dset.i))
    Dmtx.dat.ij <- Dmtx.dat[cov.dat.ij$id, cov.dat.ij$id]
    ov.ij=compute_overlap(cov.dat.ij %>% filter(Dataset == dset.i), cov.dat.ij %>% filter(Dataset == dset.j))
    
    if (length(unique(cov.dat.ij$Sex)) != 1) {
      z1 <- as.matrix(cov.dat.ij %>% dplyr::select(Sex) %>%
                        mutate(Sex=as.numeric(Sex)))
    } else {
      z1 <- NULL
    }
    if (length(unique(cov.dat.ij$Continent)) != 1) {
      z2 <- as.matrix(cov.dat.ij %>% dplyr::select(Continent) %>%
                        mutate(Continent=as.numeric(Continent)))
    } else {
      z2 <- NULL
    }
    z3 <- as.matrix(cov.dat.ij %>% dplyr::select(Age) %>%
                      mutate(Age=as.numeric(Age)))
    z <- cbind(z1, z2, z3)
    
    result <- lapply(names(fns), function(fn.name) {
      tryCatch({
        test.out <- do.call(fns[[fn.name]], list(as.dist(Dmtx.dat.ij), cov.dat.ij$Treatment, z, R=R, distance=TRUE))
        return(data.frame(Data="Raw", Method=fn.name, Dataset.Trt=dset.i,
                   Dataset.Ctrl=dset.j, Effect=test.out$statistic,
                   p.value=test.out$p.value, Overlap=ov.ij))
      }, function(e) {
        return(data.frame(Data="Raw", Method=fn.name, Dataset.Trt=dset.i,
                   Dataset.Ctrl=dset.j, Effect=NA,
                   p.value=NA, Overlap=ov.ij))
      })
    })
    
    # If both datasets are NKI, Uncorrected == Causal Crossover
    if (grepl("NKI24", dset.i) & grepl("NKI24", dset.j)) {
      result$causal <- result$uncor
      result$causal$Data <- "Crossover"
    } else {
      result$causal <- data.frame(Data="Raw", Method="Crossover", Dataset.Trt=dset.i,
                                  Dataset.Ctrl=dset.j, Effect=NA, p.value=NA,
                                  Overlap=ov.ij)
    }
    
    return(rbind(result))
  })
}, mc.cores=ncores))

saveRDS(output, "../data/pairwise.rds")