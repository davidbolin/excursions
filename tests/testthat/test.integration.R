context("Integration")
test_that("Integration L", {
  n = 11
  Q.x = Matrix(toeplitz(c(1, -0.1, rep(0, n-2))))
  L = chol(Q.x)

  seed = c(1750459768, 1598840523, 249795150, 
           2124039812, 1116456203, 1714779982)
  prob1 = gaussint(Q.chol = L, a = rep(-3,n), b = rep(3,n),
                                 seed = seed)

  expect_equal(prob1$P[1],0.9679946,tolerance=1e-7)
  expect_equal(prob1$E[1],5.986766e-06,tolerance=1e-7)
})

test_that("Integration Q", {
  n = 11
  Q.x = Matrix(toeplitz(c(1, -0.1, rep(0, n-2))))
  
  seed = c(1750459768, 1598840523, 249795150, 
           2124039812, 1116456203, 1714779982)
  prob1 = gaussint(Q = Q.x, a = rep(-3,n), b = rep(3,n),
                                 seed = seed)

  expect_equal(prob1$P[1],0.9679946,tolerance=1e-7)
  expect_equal(prob1$E[1],5.986766e-06,tolerance=1e-7)
})

test_that("Integration mu", {
  n = 11
  Q.x = Matrix(toeplitz(c(1, -0.1, rep(0, n-2))))
  mu.x = seq(-5, 5, length=n)
  seed = c(1750459768, 1598840523, 249795150, 
           2124039812, 1116456203, 1714779982)
  prob1 = gaussint(Q = Q.x, mu = mu.x, a = rep(-3,n) + mu.x, 
                                 b = rep(3,n) + mu.x, seed = seed)

  expect_equal(prob1$P[1],0.9679946,tolerance=1e-7)
  expect_equal(prob1$E[1],5.986766e-06,tolerance=1e-7)
})

test_that("Integration limit", {
  n = 11
  Q.x = Matrix(toeplitz(c(1, -0.1, rep(0, n-2))))
  mu.x = seq(-5, 5, length=n)
  seed = c(1750459768, 1598840523, 249795150, 
           2124039812, 1116456203, 1714779982)
  prob1 = gaussint(Q = Q.x, mu = mu.x, a = rep(-3,n) + mu.x, 
                                 b = rep(3,n) + mu.x, seed = seed, lim = 0.97)

  prob2 = gaussint(Q = Q.x, mu = mu.x, a = rep(-3,n) + mu.x, 
                                 b = rep(3,n) + mu.x, seed = seed, lim = 0.9)

  expect_equal(prob1$P[1],0.0,tolerance=1e-7)
  expect_equal(prob2$P[1],0.9679946,tolerance=1e-7)
})
