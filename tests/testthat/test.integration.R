context("Integration")
test_that("Integration L", {
  data <- integration.testdata1()
  prob1 <- gaussint(Q.chol = data$L, a = data$a, b = data$b, seed = data$seed,
                    max.threads = 1)
  expect_equal(prob1$P[1], 0.9679946, tolerance=10*prob1$E[1])
  #expect_equal(prob1$E[1], 5.986766e-06, tolerance=1e-6)
})

test_that("Integration Q", {
  data <- integration.testdata1()
  prob1 <- gaussint(Q = data$Q, a = data$a, b = data$b, seed = data$seed,
                    max.threads = 1)
  expect_equal(prob1$P[1], 0.9679946, tolerance=10*prob1$E[1])
  #expect_equal(prob1$E[1], 5.986766e-06, tolerance=1e-6)
})

test_that("Integration mu", {
  data <- integration.testdata1()
  prob1 <- gaussint(Q = data$Q, mu = data$mu, a = data$a + data$mu,
                    b = data$b + data$mu, seed = data$seed,
                    max.threads = 1)
  expect_equal(prob1$P[1], 0.9679946, tolerance=10*prob1$E[1])
  #expect_equal(prob1$E[1], 5.986766e-06, tolerance=1e-6)
})

test_that("Integration limit", {
  data <- integration.testdata1()
  prob1 <- gaussint(Q = data$Q, a = data$a, b = data$b,
                    seed = data$seed, lim = 0.97,
                    max.threads = 1)

  prob2 <- gaussint(Q = data$Q, a = data$a, b = data$b,
                   seed = data$seed, lim = 0.9,
                   max.threads = 1)

  expect_equal(prob1$P[1], 0.0, tolerance=10*prob1$E[1])
  #expect_equal(prob2$P[1], 0.9679946, tolerance=1e-6)
})


test_that("Integration seed", {

 seed = c(1750459768, 1598840523, 249795150,
           2124039812, 1116456203, 1714779982)
 n=10
 x <- excursions.rand(n,seed,n.threads=1)

 r <- c(0.3984955, 0.9941428, 0.3475395, 0.9096741, 0.9514461,
        0.6353268, 0.8262468, 0.7140577, 0.5236995, 0.3770943)

  expect_equal(x, r, tolerance=1e-4)
})
