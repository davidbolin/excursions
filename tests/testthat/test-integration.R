test_that("Integration L", {
  data <- integration.testdata1()
  prob1 <- gaussint(Q.chol = data$L, a = data$a, b = data$b,
                    seed = data$seed, max.threads = 1)
  expect_equal(prob1$P[1], 0.9680023, tolerance=1e-7)
  expect_equal(prob1$E[1], 5.914764e-06, tolerance=1e-6)
})

test_that("Integration Q", {
  data <- integration.testdata1()
  prob1 <- gaussint(Q = data$Q, a = data$a, b = data$b,
                    seed = data$seed, max.threads = 1)
  expect_equal(prob1$P[1], 0.9680023, tolerance=1e-7)
  expect_equal(prob1$E[1], 5.914764e-06, tolerance=1e-6)
})

test_that("Integration mu", {
  data <- integration.testdata1()
  prob1 <- gaussint(Q = data$Q, mu = data$mu, a = data$a + data$mu,
                    b = data$b + data$mu, seed = data$seed,
                    max.threads = 1)
  expect_equal(prob1$P[1], 0.9680023, tolerance=1e-7)
  expect_equal(prob1$E[1], 5.914764e-06, tolerance=1e-6)
})

test_that("Integration limit", {
  data <- integration.testdata1()
  prob1 <- gaussint(Q = data$Q, a = data$a, b = data$b,
                    seed = data$seed, lim = 0.97,
                    max.threads = 1)

  prob2 <- gaussint(Q = data$Q, a = data$a, b = data$b,
                   seed = data$seed, lim = 0.9,
                   max.threads = 1)

  expect_equal(prob1$P[1], 0.0, tolerance=1e-7)
  expect_equal(prob2$P[1],  0.9680023, tolerance=1e-6)
})

test_that("Integration reordering", {
  data <- integration.testdata1()
  prob1 <- gaussint(Q = data$Q, mu = data$mu, a = data$a + data$mu,
                    b = data$b + data$mu, seed = data$seed,
                    max.threads = 1,use.reordering = "sparsity")
  expect_equal(prob1$P[1], 0.9680023, tolerance=1e-5)
  expect_equal(prob1$E[1], 5.914764e-06, tolerance=1e-5)
})
