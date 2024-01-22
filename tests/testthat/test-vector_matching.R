test_that("one-hot encoding works with numeric vector", {
  Ts <- c(1, 1, 0, 0)
  Ts.ohe <- as.matrix(ohe(Ts))
  expect_true(all(Ts.ohe[1,] == Ts.ohe[2,]))
  expect_true(all(Ts.ohe[3,] == Ts.ohe[4,]))
  expect_false(all(Ts.ohe[1,] == Ts.ohe[3,]))
})

test_that("one-hot encoding works with character vector", {
  Ts <- c("a", "a", "b", "b")
  Ts.ohe <- as.matrix(ohe(Ts))
  expect_true(all(Ts.ohe[1,] == Ts.ohe[2,]))
  expect_true(all(Ts.ohe[3,] == Ts.ohe[4,]))
  expect_false(all(Ts.ohe[1,] == Ts.ohe[3,]))
})

test_that("one-hot encoding works with string vector", {
  Ts <- c("aaa", "aaa", "bbb", "bbb")
  Ts.ohe <- as.matrix(ohe(Ts))
  expect_true(all(Ts.ohe[1,] == Ts.ohe[2,]))
  expect_true(all(Ts.ohe[3,] == Ts.ohe[4,]))
  expect_false(all(Ts.ohe[1,] == Ts.ohe[3,]))
})

test_that("one-hot encoding works with factor vector", {
  Ts <- factor(c("aaa", "aaa", "bbb", "bbb"))
  Ts.ohe <- as.matrix(ohe(Ts))
  expect_true(all(Ts.ohe[1,] == Ts.ohe[2,]))
  expect_true(all(Ts.ohe[3,] == Ts.ohe[4,]))
  expect_false(all(Ts.ohe[1,] == Ts.ohe[3,]))
})

test_that("zero-one distance computation works", {
  Ts <- c(1, 1, 0, 0)
  DT <- as.matrix(zero_one_dist(Ts))
  expect_true(all(DT[1:2, 1:2] == 0))
  expect_true(all(DT[3:4, 3:4] == 0))
  expect_true(all(DT[1:2, 3:4] == 1))
  expect_true(all(DT[3:4, 1:2] == 1))
})

test_that("zero-one distance computation works with arbitrary ordering and string", {
  Ts <- c("aaa", "bbb", "aaa", "bbb", "ccc")
  DT <- as.matrix(zero_one_dist(Ts))
  DT.true <- cbind(c(0, 1, 0, 1, 1), c(1, 0, 1, 0, 1), c(0, 1, 0, 1, 1),
                   c(1, 0, 1, 0, 1), c(1, 1, 1, 1, 0))
  expect_true(all(DT == DT.true))
})

test_that("vector matching with one odd-ball per-group", {
  nrep <- 100
  
  res <- sapply(1:nrep, function(i) {
    Ts <- c(rep(1, 100), rep(2, 100))
    Xs <- cbind(c(-4, runif(99), runif(99), 5), runif(200), runif(200))
    
    retained.ids <- suppressWarnings(cb.align.vm_trim(Ts, Xs))
    
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

test_that("moderate overlap example retains more samples than limited overlap example", {
  
})