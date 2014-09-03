context("Integration")
test_that("Integration", {
library(Matrix)

n = 11
mu.x = seq(-5, 5, length=n)
Q.x = Matrix(toeplitz(c(1, -0.1, rep(0, n-2))))
L = chol(Q.x)

seed = c(1750459768, 1598840523, 249795150, 2124039812, 1116456203, 1714779982)
prob1 = excursions.integration(L, a = rep(-3,n), b = rep(3,n),seed = seed)

expect_equal(prob1$P[1],0.9679946,tolerance=1e-7)
expect_equal(prob1$E[1],5.986766e-06,tolerance=1e-7)
})
