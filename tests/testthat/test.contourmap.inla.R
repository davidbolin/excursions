context("Contourmap.inla")

test_that("Contourmap.inla, stack extraction", {
  if (require.nowarnings("INLA")) {
    data <- testdata.inla()

    res1 = contourmap.inla(data$result, data$stack, tag = "pred",
                           n.levels=4,seed=data$seed,
                           compute = list(F = FALSE))

    mu <- c(-2.0987990,-0.4124901,1.6699354,2.6138399,4.8317157,4.8170408,
            1.7105803,-0.4731774,0.8455982,-0.6867296,-0.4645699)
    expect_equal(res1$meta$mu,mu,tolerance=1e-7)
  }
})


test_that("Contourmap.inla, test ind", {
  if (require.nowarnings("INLA")) {
  data <- testdata.inla()
  ind1 = c(1,2,3,4)
  ind2 = c(4,3,2,1)
  ind3 = rep(FALSE,data$n)
  ind3[1:4] = TRUE

  res1 = contourmap.inla(data$result, data$stack, tag = "pred",
                         n.levels=2,ind=ind1, seed=data$seed,alpha=0.1,
                         max.threads=1)
  res2 = contourmap.inla(data$result, data$stack, tag = "pred",
                         n.levels=2,ind=ind2, seed=data$seed,alpha=0.1,
                         max.threads=1)
  res3 = contourmap.inla(data$result, data$stack, tag = "pred",
                         n.levels=2,ind=ind3, seed=data$seed,alpha=0.1,
                         max.threads=1)

  expect_equal(res1$F,res2$F,tolerance=1e-7)
  expect_equal(res2$F,res3$F,tolerance=1e-7)
}
})

test_that("Contourmap.inla, P measures", {
  if (require.nowarnings("INLA")) {
    data <- testdata.inla()

    res1 = contourmap.inla(data$result, data$stack, tag = "pred",
                           n.levels=4,seed=data$seed,
                           max.threads=1,
                           compute = list(F = FALSE, measures = c("P2","P1")))

    expect_equal(res1$P1,0.9760582,tolerance=1e-7)
    expect_equal(res1$P2,0.6282216,tolerance=1e-7)
  }
})
