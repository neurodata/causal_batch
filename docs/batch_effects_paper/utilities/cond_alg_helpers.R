require(energy)
require(cdcsis)
require(causalBatch)
require(sva)
require(GeneralisedCovarianceMeasure)

# Conditional distance correlation from Wang et al., 2015
cond.dcorr <- function(Ys, Ts, Xs, R=1000, dist.method="euclidean", distance = FALSE, seed=1, num.threads=1) {
  Xs <- as.data.frame(Xs)
  
  if (isTRUE(distance)) {
    DY <- as.dist(Ys)
  } else {
    DY = dist(Ys, method=dist.method)
  }
  
  DT = causalBatch:::zero_one_dist(Ts)
  
  # run statistical test
  test.out <- cdcov.test(DY, DT, Xs, num.bootstrap = R,
                         seed=seed, num.threads=num.threads, distance=TRUE)
  return(list(Test=test.out))
}

dcorr <- function(Ys, Ts, Xs, R=1000, dist.method="euclidean", distance = FALSE, seed=1, num.threads=1) {
  Xs <- as.data.frame(Xs)
  
  if (isTRUE(distance)) {
    DY <- as.dist(Ys)
  } else {
    DY = dist(Ys, method=dist.method)
  }
  
  DT = causalBatch:::zero_one_dist(Ts)
  
  test.out <- dcor.test(DY, DT, R=R)
  return(list(Test=test.out))
}

cond.combat <- function(Ys, Ts, Xs, ...) {
  mod <- model.matrix(as.formula("~factor(Sex) + Age"), data=Xs)
  Ys.cor <- t(ComBat(t(Ys), Ts, mod=mod))
  return(list(Ys.corrected=Ys.cor, Ts=Ts, Xs=Xs))
}

assoc.combat <- function(Ys, Ts, Xs, match.form, match.args=NULL, retain.ratio=0.05, ...) {
  Ys.cor <- t(ComBat(t(Ys), Ts))
  return(list(Ys.corrected=Ys.cor, Ts=Ts, Xs=Xs))
}

raw.preproc <- function(Ys, Ts, Xs, match.form, match.args=NULL, retain.ratio=0.05, ...) {
  return(list(Ys.corrected=Ys, Ts=Ts, Xs=Xs))
}


gcm <- function(Y, G, X, R=1000, regr.method="xgboost") {
  res <- tryCatch({
    gcm.test(as.matrix(Y), as.matrix(G), as.matrix(X), regr.method=regr.method, nsim=R)},
    error=function(e) {
      return(list(stat=NaN, pvalue=NaN))
    })
  return(list(statistic=res$test.statistic, p.value=res$p.value))
}

matching.combat <- function(Ys, Ts, Xs, match.form, match.args=NULL, retain.ratio=0.05, ...) {
  return(cb.correct.matching_cComBat(Ys, Ts, Xs, match.form, match.args=match.args, retain.ratio=retain.ratio))
}

compute_overlap <- function(X1, X2) {
  # probability of drawing two individuals with the same sex
  X1.sex <- X1 %>%
    group_by(Sex) %>%
    summarize(Per=n(), .groups="keep") %>%
    ungroup() %>%
    mutate(Per=Per/sum(Per))
  X2.sex <- X2 %>%
    group_by(Sex) %>%
    summarize(Per=n(), .groups="keep") %>%
    ungroup() %>%
    mutate(Per=Per/sum(Per))
  range.age <- c(min(X1$Age, X2$Age), max(X1$Age, X2$Age))
  per.ov <- sum(sapply(1:2, function(sex) {
    tryCatch({
      ov.sex.X1 <- (X1.sex %>% filter(Sex == sex))$Per
      ov.sex.X2 <- (X2.sex %>% filter(Sex == sex))$Per
      X1.sex.age <- (X1 %>% filter(Sex == sex) %>% select(Age))$Age
      X2.sex.age <- (X2 %>% filter(Sex == sex) %>% select(Age))$Age
      # obtain pdf for age
      X1.dens <- density(as.numeric(X1.sex.age), from=min(range.age), to=max(range.age))
      X1.dens$y <- X1.dens$y/sum(X1.dens$y)
      X2.dens <- density(as.numeric(X2.sex.age), from=min(range.age), to=max(range.age))
      X2.dens$y <- X2.dens$y/sum(X2.dens$y)
      ov.sex.age <- sum(pmin(X1.dens$y*ov.sex.X1, X2.dens$y*ov.sex.X2))
      return(ov.sex.age)
    }, error=function(e) {return(0)})
  }))
  return(as.numeric(unique((X1$Continent)) == unique(X2$Continent))*per.ov)
}
