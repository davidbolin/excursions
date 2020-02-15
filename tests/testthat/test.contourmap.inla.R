context("INLA interface")

test_that("stack extraction", {
  if (requireNamespace("INLA", quietly = TRUE)) {
    data <- testdata.inla.small()
    ind <- excursions:::inla.output.indices(data$result,stack=data$stack,tag="pred")
    expect_equal(ind,c(6,7,8,9,10,11,12),tolerance=1e-7)
    for(i in 1:data$result$misc$configs$nconfig){
      config = excursions:::private.get.config(data$result,i)
      if(config$lp == 0)
        break
    }
    mu <- c(3.2283547,3.2631704,1.7001027,-0.4210484,1.0174016,-1.4767174,-0.4210482,
            3.2283544,3.2631702 ,1.0174017,1.7001023,-1.8092808,-0.4210482,3.2283544,
            3.2631702,1.0174017,1.7001023,-1.4767174,-1.8092808,-0.2789614,-0.8202960,
            1.2320903,1.2669076,0.6181522,2.8978491,0.9854832,1.1179056)
    expect_equal(config$mu,mu,tolerance=1e-2)
    vars <- c(5.028681e-05,5.028681e-05,5.028692e-05,5.028659e-05,5.028642e-05,2.255220e+00,
              5.089838e-05,5.089860e-05,5.089860e-05,5.089821e-05,5.089872e-05,3.240656e+00,
              5.059248e-05,5.059269e-05,5.059269e-05,5.059231e-05,5.059281e-05,2.255219e+00,
              3.240656e+00,1.812985e+00,4.011360e-02,1.001390e+00,1.001390e+00,4.011075e-02,
              3.605275e-01,1.850041e+00,3.140320e-01)
    expect_equal(config$vars,vars,tolerance=1e-2)
  }
})

context("Contourmap.inla")

test_that("Contourmap.inla, test ind", {
  if (requireNamespace("INLA", quietly = TRUE)) {
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

  expect_equal(res1$F,res2$F,tolerance=1e-4)
  expect_equal(res2$F,res3$F,tolerance=1e-4)
}
})

test_that("Contourmap.inla, P measures", {
  if (requireNamespace("INLA", quietly = TRUE)) {
    data <- testdata.inla.small()

    res1 = contourmap.inla(data$result, data$stack, tag = "pred",
                           n.levels=4,seed=data$seed,
                           max.threads=1,
                           compute = list(F = FALSE, measures = c("P2","P1")),
                           method='EB')

    expect_equal(res1$P1,0.7732031,tolerance=1e-3)
    expect_equal(res1$P2,0.6558154,tolerance=1e-3)

    res1 = contourmap.inla(data$result, data$stack, tag = "pred",
                           n.levels=4,seed=data$seed,
                           max.threads=1,
                           compute = list(F = FALSE, measures = c("P2","P1")),
                           method='QC')
    expect_equal(res1$P1,0.7669963,tolerance=1e-3)
    expect_equal(res1$P2,0.6613112,tolerance=1e-3)
  }
})
