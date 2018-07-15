context("Contourmap.inla")

test_that("Contourmap.inla, stack extraction", {
  if (require.nowarnings("INLA")) {
    data <- testdata.inla.small()

    res1 = contourmap.inla(data$result, data$stack, tag = "pred",
                           n.levels=4,seed=data$seed,
                           compute = list(F = FALSE))

    mu <- c(-1.3290310,-0.4210768,3.2283824,3.2631843,1.0173736,1.7001533,-3.7534800)
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
    data <- testdata.inla.small()

    res1 = contourmap.inla(data$result, data$stack, tag = "pred",
                           n.levels=4,seed=data$seed,
                           max.threads=1,
                           compute = list(F = FALSE, measures = c("P2","P1")))

    expect_equal(res1$P1,0.7395756,tolerance=1e-7)
    expect_equal(res1$P2,0.4400722,tolerance=1e-7)
  }
})
