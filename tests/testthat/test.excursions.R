context("Excursions")


test_that("Excursions, alpha = 1, type = >", {
  data <- integration.testdata1()
  res <- excursions(alpha=1, u=0, mu=data$mu, Q=data$Q, type='>',
                    seed = data$seed,
                    max.threads = 1)
  r <- c(2.467173e-15, 1.031418e-09, 7.755949e-06, 2.543606e-03,
       7.599777e-02, 4.194159e-01, 8.192652e-01, 9.747060e-01,
       9.984658e-01, 9.999621e-01, 9.999997e-01)
  expect_equal(res$F,r,tolerance=1e-7)
})

test_that("Excursions, alpha = 1, type = <", {
  data <- integration.testdata1()
  res = excursions(alpha=1, u=0, mu=data$mu, Q=data$Q, type='<', seed = data$seed,
                    max.threads = 1)
  r = c(9.999997e-01, 9.999622e-01, 9.984766e-01, 9.746239e-01, 8.190389e-01,
       4.192016e-01, 7.586059e-02, 2.540068e-03, 7.753973e-06, 1.033270e-09,
       2.468060e-15)
  expect_equal(res$F,r,tolerance=1e-7)
})

test_that("Excursions, alpha = 1, type = =", {
  data <- integration.testdata1()
  res = excursions(alpha=1, u=0, mu=data$mu+0.1, Q=data$Q, type="=", seed = data$seed,
                    max.threads = 1)
  r = c(7.381175e-07, 8.177720e-05, 3.204621e-03, 5.127762e-02, 3.327196e-01,
       6.420550e-01, 1.812984e-01, 2.193460e-02, 1.155138e-03, 2.547679e-05,
       1.945911e-07)
  expect_equal(res$F,r,tolerance=1e-7)
})

test_that("Excursions, alpha = 1, type = !=", {
  data <- integration.testdata1()
  res = excursions(alpha=1, u=0, mu=data$mu+0.1, Q=data$Q, type='!=', seed = data$seed,
                    max.threads = 1)
  r = c(0.9999993, 0.9999182, 0.9967954, 0.9487224, 0.6672804, 0.3579450,
       0.8187016, 0.9780654, 0.9988449, 0.9999745, 0.9999998)
  expect_equal(res$F,r,tolerance=1e-7)
})

test_that("Excursions, alpha = 0.1, type = >", {
  data <- integration.testdata1()
  res = excursions(alpha=0.1, u=0, mu=data$mu+0.1, Q=data$Q, type='>', seed = data$seed,
                    max.threads = 1)
  r = c(0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000,
       0.0000000, 0.9801446, 0.9988957, 0.9999751, 0.9999998)
  res$F[is.na(res$F)] = 0
  expect_equal(res$F,r,tolerance=1e-7)
})

test_that("Excursions, alpha = 0.1, type = <", {
  data <- integration.testdata1()
  res = excursions(alpha=0.1, u=0, mu=data$mu+0.1, Q=data$Q, type='<', seed = data$seed,
                    max.threads = 1)
  r = c(0.9999429, 0.9999434, 0.9978976, 0.9679410, 0.0000000, 0.0000000,
       0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000)
  res$F[is.na(res$F)] = 0
  expect_equal(res$F,r,tolerance=1e-7)
})

test_that("Excursions, alpha = 0.1, type = =", {
  data <- integration.testdata1()
  res = excursions(alpha=0.1, u=0, mu=data$mu+0.1, Q=data$Q, type='=', seed = data$seed,
                    max.threads = 1)
  r = c(7.381175e-07, 8.177720e-05, 3.204621e-03, 5.127762e-02, 1.000000e+00,
       1.000000e+00, 1.000000e+00, 2.193460e-02, 1.155138e-03, 2.547679e-05,
       1.945911e-07)
  res$F[is.na(res$F)] = 1
  expect_equal(res$F,r,tolerance=1e-7)
})

test_that("Excursions, alpha = 0.1, type = !=", {
  data <- integration.testdata1()
  res = excursions(alpha=0.1, u=0, mu=data$mu+0.1, Q=data$Q, type='!=', seed = data$seed,
                    max.threads = 1)
  r = c(0.9999993, 0.9999182, 0.9967954, 0.9487224, 0.0000000, 0.0000000,
       0.0000000, 0.9780654, 0.9988449, 0.9999745, 0.9999998)
  res$F[is.na(res$F)] = 0
  expect_equal(res$F,r,tolerance=1e-7)
})

test_that("Excursions, move u to mu", {
  data <- integration.testdata1()

  res = excursions(alpha=0.1, u=1, mu=data$mu, Q=data$Q, type='>', seed = data$seed,
                    max.threads = 1)
  res2 = excursions(alpha=0.1, u=0, mu=data$mu-1, Q=data$Q, type='>', seed = data$seed,
                    max.threads = 1)
  res$F[is.na(res$F)] = 0
  res2$F[is.na(res2$F)] = 0
  expect_equal(res$F,res2$F,tolerance=1e-7)
})

test_that("Excursions, input variances", {
  data <- integration.testdata1()

  vars = diag(solve(data$Q))
  res1 = excursions(alpha=1, u=0, mu=data$mu, Q=data$Q, type='>', seed = data$seed, vars = vars,
                    max.threads = 1)
  res2 = excursions(alpha=1, u=0, mu=data$mu, Q=data$Q, type='>',seed = data$seed,
                    max.threads = 1)
  expect_equal(res1$F,res2$F,tolerance=1e-7)
})

test_that("Excursions, ind argument order", {
  data <- integration.testdata1()

  vars = diag(solve(data$Q))

  ind1 = c(1,2,3,4)
  ind2 = c(4,3,2,1)
  ind3 = rep(FALSE,length(data$mu))
  ind3[1:4] = TRUE
  res1 = excursions(alpha=1, u=0, mu=data$mu, Q=data$Q, type='>', seed = data$seed, ind = ind1,
                    max.threads = 1)
  res2 = excursions(alpha=1, u=0, mu=data$mu, Q=data$Q, type='>', seed = data$seed, ind = ind2,
                    max.threads = 1)
  res3 = excursions(alpha=1, u=0, mu=data$mu, Q=data$Q, type='>', seed = data$seed, ind = ind3,
                    max.threads = 1)

  expect_equal(res1$F,res2$F,tolerance=1e-7)
  expect_equal(res2$F,res3$F,tolerance=1e-7)
})

#Tests to add:

#Test that Q.chol and Q gives the same result

#Test the max.size argument

#Test reo

#Test rho

#Test max.threads

#Test QC method

#Test ind argument


#res1 = excursions2::excursions(alpha=0.1, u=1, mu=mu.x+0.1, Q=Q.x, type='!=', seed = seed, max.threads = 1)
#res2 = excursions::excursions(alpha=0.1, u=1, mu=mu.x+0.1, Q=Q.x, type='!=', max.threads = 1)

#plot(res1$F)
#lines(res2$F, col=2)
