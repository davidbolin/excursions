test_that("stack extraction", {
  skip_on_cran()
  local_exc_safe_inla()

    data <- testdata.inla.small()
    ind <- excursions:::inla.output.indices(data$result,stack=data$stack,tag="pred")
    expect_equal(ind,c(6,7,8,9,10,11,12),tolerance=1e-7)
    for(i in 1:data$result$misc$configs$nconfig){
      config = excursions:::private.get.config(data$result,i)
      if(config$lp == 0)
        break
    }

    # Only check prediction of unobserved values
    expect_snapshot_value(config$mu[ind[2:6]],
                          style = "serialize",
                          tolerance = 1e-2)

})

test_that("Contourmap.inla, test ind", {
  skip_on_cran()
  local_exc_safe_inla()

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

})

test_that("Contourmap.inla, P measures", {
  skip_on_cran()
  local_exc_safe_inla()

    data <- testdata.inla.small()

    ind <- 2:6
    res1 = contourmap.inla(data$result, data$stack, tag = "pred",
                           n.levels=4,seed=data$seed,
                           max.threads=1,
                           ind = ind,
                           compute = list(F = FALSE, measures = c("P2","P1")),
                           method='EB')

    expect_equal(res1$P1,0.9963946,tolerance=2e-2)
    expect_equal(res1$P2,0.8983837,tolerance=2e-2)

    res2 = contourmap.inla(data$result, data$stack, tag = "pred",
                           n.levels=4,seed=data$seed,
                           max.threads=1,
                           ind = ind,
                           compute = list(F = FALSE, measures = c("P2","P1")),
                           method='QC')
    expect_equal(res2$P1,0.9953097,tolerance=2e-2)
    expect_equal(res2$P2,0.8893297,tolerance=2e-2)

})
