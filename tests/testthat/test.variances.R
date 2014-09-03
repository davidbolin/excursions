context("Variances")
test_that("Variances", {
library(Matrix)

n = 11
mu.x = seq(-5, 5, length=n)
Q.x = Matrix(toeplitz(c(1, -0.1, rep(0, n-2))))
L = chol(Q.x)

vars <- excursions.variances(L)
v = diag(solve(Q.x))
expect_equal(vars,v,tolerance=1e-7)
})