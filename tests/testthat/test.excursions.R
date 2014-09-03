context("Excursions")
test_that("Excursions", {
n = 11
Q.x = sparseMatrix(i=c(1:n, 2:n),
                   j=c(1:n, 1:(n-1)),
                   x=c(rep(1, n), rep(-0.1, n-1)),
                   dims=c(n, n),
                   symmetric=TRUE)

mu.x = seq(-5, 5, length=n)

seed = c(1750459768, 1598840523, 249795150, 2124039812, 1116456203, 1714779982)

# Tests with alpha = 1
res1 = excursions(alpha=1, u=0, mu=mu.x, Q=Q.x, type='>', seed = seed)
r1 = c(2.467173e-15, 1.031418e-09, 7.755949e-06, 2.543606e-03, 
       7.599777e-02, 4.194159e-01, 8.192652e-01, 9.747060e-01, 
       9.984658e-01, 9.999621e-01, 9.999997e-01)
expect_equal(res1$F,r1,tolerance=1e-7)

res2 = excursions(alpha=1, u=0, mu=mu.x, Q=Q.x, type='<', seed = seed)
r2 = c(9.999997e-01, 9.999622e-01, 9.984766e-01, 9.746239e-01, 8.190389e-01,
       4.192016e-01, 7.586059e-02, 2.540068e-03, 7.753973e-06, 1.033270e-09,     
       2.468060e-15)
expect_equal(res2$F,r2,tolerance=1e-7)

res3 = excursions(alpha=1, u=0, mu=mu.x+0.1, Q=Q.x, type="=", seed = seed)
r3 = c(7.381175e-07, 8.177720e-05, 3.204621e-03, 5.127762e-02, 3.327196e-01,
       6.420550e-01, 1.812984e-01, 2.193460e-02, 1.155138e-03, 2.547679e-05,
       1.945911e-07)
expect_equal(res3$F,r3,tolerance=1e-7)

res4 = excursions(alpha=1, u=0, mu=mu.x+0.1, Q=Q.x, type='!=', seed = seed)
r4 = c(0.9999993, 0.9999182, 0.9967954, 0.9487224, 0.6672804, 0.3579450,
       0.8187016, 0.9780654, 0.9988449, 0.9999745, 0.9999998)
expect_equal(res4$F,r4,tolerance=1e-7)

# Tests with alpha = 0.1
res5 = excursions(alpha=0.1, u=0, mu=mu.x+0.1, Q=Q.x, type='>', seed = seed)
r5 = c(0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 
       0.0000000, 0.9801446, 0.9988957, 0.9999751, 0.9999998)
expect_equal(res5$F,r5,tolerance=1e-7)

res6 = excursions(alpha=0.1, u=0, mu=mu.x+0.1, Q=Q.x, type='<', seed = seed)
r6 = c(0.9999429, 0.9999434, 0.9978976, 0.9679410, 0.0000000, 0.0000000,
       0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000)
expect_equal(res6$F,r6,tolerance=1e-7)

res7 = excursions(alpha=0.1, u=0, mu=mu.x+0.1, Q=Q.x, type='=', seed = seed)
r7 = c(7.381175e-07, 8.177720e-05, 3.204621e-03, 5.127762e-02, 1.000000e+00, 
       1.000000e+00, 1.000000e+00, 2.193460e-02, 1.155138e-03, 2.547679e-05,
       1.945911e-07)
expect_equal(res7$F,r7,tolerance=1e-7)

res8 = excursions(alpha=0.1, u=0, mu=mu.x+0.1, Q=Q.x, type='!=', seed = seed) 
r8 = c(0.9999993, 0.9999182, 0.9967954, 0.9487224, 0.0000000, 0.0000000, 
       0.0000000, 0.9780654, 0.9988449, 0.9999745, 0.9999998)
expect_equal(res8$F,r8,tolerance=1e-7)

#Test that u can be moved to mu
res9 = excursions(alpha=0.1, u=1, mu=mu.x, Q=Q.x, type='>', seed = seed)
res10 = excursions(alpha=0.1, u=0, mu=mu.x-1, Q=Q.x, type='>', seed = seed)
expect_equal(res9$F,res10$F,tolerance=1e-7)

#Test that variances can be added
vars = diag(solve(Q.x))
res11 = excursions(alpha=1, u=0, mu=mu.x, Q=Q.x, type='>', seed = seed, vars = vars)
res12 = excursions(alpha=1, u=0, mu=mu.x, Q=Q.x, type='>',seed = seed)
expect_equal(res11$F,res12$F,tolerance=1e-7)

#Test that Q.chol and Q gives the same result

#Test the max.size argument

#Test reo

#Test rho

#Test max.threads

#Test QC method

#Test ind argument

})
 

#res1 = excursions2::excursions(alpha=0.1, u=1, mu=mu.x+0.1, Q=Q.x, type='!=', seed = seed)
#res2 = excursions::excursions(alpha=0.1, u=1, mu=mu.x+0.1, Q=Q.x, type='!=')

#plot(res1$F)
#lines(res2$F, col=2)