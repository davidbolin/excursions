context("Variances")
test_that("Variances", {
  data <- integration.testdata1()
  vars <- excursions.variances(data$L)
  v = diag(solve(data$Q))
  expect_equal(vars,v,tolerance=1e-7)
})


test_that("Variances Q and L", {
  data <- integration.testdata1()
  v1 <- excursions.variances(L = data$L)
  v2 <- excursions.variances(Q = data$Q)
  expect_equal(v1,v2,tolerance=1e-7)
})