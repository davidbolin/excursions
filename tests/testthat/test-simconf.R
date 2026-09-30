test_that("Simconf", {
  data <- integration.testdata1()
  res <- simconf(Q = data$Q, mu = data$mu, seed = data$seed, alpha = 0.1, max.threads = 1)
  ra <- c(
    -7.6070830, -6.6203520, -5.6204871, -4.6204885, -3.6204885,
    -2.6204885, -1.6204885, -0.6204885, 0.3795129, 1.3796480, 2.3929170
  )
  rb <- c(
    -2.3929170, -1.3796480, -0.3795129, 0.6204885, 1.6204885, 2.6204885,
    3.6204885, 4.6204885, 5.6204871, 6.6203520, 7.6070830
  )
  expect_equal(res$a, ra, tolerance = 1e-3)
  expect_equal(res$b, rb, tolerance = 1e-3)
})


test_that("Simultaneous confidence bands are unchanged", {
  d <- testdata.spde(6)
  r <- simconf(alpha = 0.1, mu = d$mu, Q = d$Q, seed = d$seed, max.threads = 1)
  expect_equal(
    c(a = ref.summary(r$a), b = ref.summary(r$b)),
    REF$sc[c("a.sum", "a.wsum", "b.sum", "b.wsum")],
    tolerance = 1e-8
  )
})

test_that("Simconf gives the marginal bounds in the right order", {
  ## a.marginal and b.marginal used to be swapped
  d <- testdata.spde(6)
  r <- simconf(alpha = 0.1, mu = d$mu, Q = d$Q, seed = d$seed, max.threads = 1, n.iter = 1000)
  sd <- sqrt(excursions.variances(Q = d$Q))
  expect_equal(r$a.marginal, d$mu + qnorm(0.05) * sd, tolerance = 1e-10)
  expect_equal(r$b.marginal, d$mu + qnorm(0.95) * sd, tolerance = 1e-10)
  expect_true(all(r$a < r$a.marginal & r$a.marginal < r$mean))
  expect_true(all(r$b > r$b.marginal & r$b.marginal > r$mean))
})

test_that("Simconf uses n.iter", {
  d <- testdata.spde(6)
  seen <- NULL
  gi <- gaussint
  local_mocked_bindings(gaussint = function(..., n.iter) {
    seen <<- c(seen, n.iter)
    gi(..., n.iter = n.iter)
  })
  r <- simconf(alpha = 0.1, mu = d$mu, Q = d$Q, seed = d$seed, max.threads = 1, n.iter = 500)
  expect_true(length(seen) > 0)
  expect_true(all(seen == 500))
})

test_that("Simconf with integer and logical indices agree", {
  d <- testdata.spde(6)
  ind <- 10:20
  lind <- seq_len(d$n) %in% ind
  r.int <- simconf(alpha = 0.1, mu = d$mu, Q = d$Q, ind = ind, seed = d$seed, max.threads = 1, n.iter = 2000)
  r.log <- simconf(alpha = 0.1, mu = d$mu, Q = d$Q, ind = lind, seed = d$seed, max.threads = 1, n.iter = 2000)
  r.all <- simconf(alpha = 0.1, mu = d$mu, Q = d$Q, seed = d$seed, max.threads = 1, n.iter = 2000)
  expect_equal(r.int$a, r.log$a)
  expect_equal(r.int$b, r.log$b)
  expect_length(r.int$a, length(ind))
  ## The band for fewer nodes is narrower than the band for all nodes
  expect_true(all(r.int$b - r.int$a < r.all$b[ind] - r.all$a[ind]))
  ## and wider than the marginal band
  expect_true(all(r.int$b - r.int$a > r.int$b.marginal - r.int$a.marginal))
})
