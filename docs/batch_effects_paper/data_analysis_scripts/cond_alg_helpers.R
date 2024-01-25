# Conditional distance correlation from Wang et al., 2015
cond.dcorr <- function(Ys, Ts, Xs, R=1000, dist.method="euclidean", distance = FALSE, seed=1, num.threads=1) {
  Xs <- as.data.frame(Xs)
  
  if (isTRUE(distance)) {
    DY <- as.dist(Ys)
  } else {
    DY = dist(Ys, method=dist.method)
  }
  
  DT = zero_one_dist(Ts)
  
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
  
  DT = zero_one_dist(Ts)
  
  test.out <- dcor.test(DY, DT, R=R)
  return(list(Test=test.out))
}