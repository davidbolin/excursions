context("Variances")
test_that("Variances", {
  data <- integration.testdata1()
  vars <- excursions.variances(data$L)
  v = diag(solve(data$Q)) 
  expect_equal(vars,v,tolerance=1e-7)
})