test_that("stack extraction", {
  skip_on_cran()
  local_exc_safe_inla()

    data <- testdata.inla.small()
    tmp <- excursions:::inla.output.indices(data$result, stack = data$stack,
                                            tag="pred")
    ind <- tmp$index
    if(tmp$result.updated){
        result <- tmp$result    
    } else {
        result <- data$result
    }
    
    expect_equal(ind, c(6,7,8,9,10,11,12), tolerance = 1e-7)
    for(i in 1:result$misc$configs$nconfig){
      config = excursions:::private.get.config(result, i)
      if(config$lp == 0)
        break
    }

    # Only check prediction of unobserved values
    expect_snapshot_value(config$mu[ind[2:6]],
                          style = "serialize",
                          tolerance = 1e-2)
    
    #test compact mode
    data2 <- testdata.inla.small(inla.mode = "compact")
    tmp <- excursions:::inla.output.indices(data2$result, stack = data$stack,
                                            tag="pred", compressed = FALSE)
    ind2 <- tmp$index
    result <- tmp$result    
    expect_equal(ind2, c(6,7,8,9,10,11,12), tolerance = 1e-7)
    for(i in 1:result$misc$configs$nconfig){
        config2 = excursions:::private.get.config(result,i)
        if(config2$lp == 0)
            break
    }
    expect_equal(config2$mu, config$mu, tolerance = 1e-2)

    tmp <- excursions:::inla.output.indices(data2$result, stack = data$stack,
                                            tag="pred", compressed = TRUE)
    ind3 <- tmp$index
    result <- tmp$result    
    expect_equal(ind3, 1:length(c(6,7,8,9,10,11,12)), tolerance = 1e-7)
    for(i in 1:result$misc$configs$nconfig){
        config3 = excursions:::private.get.config(result,i)
        if(config3$lp == 0)
            break
    }
    expect_equal(config3$mu[ind3], config2$mu[ind2], tolerance = 1e-2)
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
                         n.levels = 2, ind = ind1, seed=data$seed, 
                         alpha = 0.1, max.threads = 1)
  res2 = contourmap.inla(data$result, data$stack, tag = "pred",
                         n.levels = 2, ind = ind2, seed = data$seed,
                         alpha = 0.1, max.threads = 1)
  res3 = contourmap.inla(data$result, data$stack, tag = "pred",
                         n.levels = 2, ind = ind3, seed = data$seed,
                         alpha = 0.1, max.threads = 1)

  expect_equal(res1$F, res2$F, tolerance = 1e-4)
  expect_equal(res2$F, res3$F, tolerance = 1e-4)

  data2 <- testdata.inla(inla.mode = "compact")
  res4 = contourmap.inla(data2$result, data2$stack, tag = "pred",
                         n.levels = 2, ind = ind1, seed = data2$seed,
                         alpha = 0.1, max.threads = 1, compressed = FALSE)
  res5 = contourmap.inla(data2$result, data2$stack, tag = "pred",
                         n.levels = 2, ind = ind2, seed=data2$seed,
                         alpha = 0.1, max.threads = 1, compressed = FALSE)
  res6 = contourmap.inla(data2$result, data2$stack, tag = "pred",
                         n.levels = 2, ind = ind3, seed = data2$seed,
                         alpha = 0.1, max.threads = 1, compressed = FALSE)
  
  expect_equal(res3$F[1:4], res4$F[1:4], tolerance = 5e-2)
  expect_equal(res4$F, res5$F, tolerance = 1e-4)
  expect_equal(res5$F, res6$F, tolerance = 1e-4)
  
  res7 = contourmap.inla(data2$result, data2$stack, tag = "pred",
                         n.levels = 2, ind = ind1, seed = data2$seed,
                         alpha = 0.1, max.threads = 1, compressed = TRUE)
  res8 = contourmap.inla(data2$result, data2$stack, tag = "pred",
                         n.levels = 2, ind = ind2, seed=data2$seed,
                         alpha = 0.1, max.threads = 1, compressed = TRUE)
  res9 = contourmap.inla(data2$result, data2$stack, tag = "pred",
                         n.levels = 2, ind = ind3, seed = data2$seed,
                         alpha = 0.1, max.threads = 1, compressed = TRUE)
  
  expect_equal(res6$F[1:4], res7$F[1:4], tolerance = 5e-2)
  expect_equal(res7$F, res8$F, tolerance = 1e-4)
  expect_equal(res8$F, res9$F, tolerance = 1e-4)
  
})


test_that("Contourmap.inla, P measures", {
  skip_on_cran()
  local_exc_safe_inla()

    data <- testdata.inla.small()

    ind <- 2:6
    res1 = contourmap.inla(data$result, data$stack, tag = "pred",
                           n.levels = 4, seed = data$seed,
                           max.threads = 1, ind = ind,
                           compute = list(F = FALSE, measures = c("P2","P1")),
                           method='EB')

    expect_snapshot_value(res1$P1,
                          style = "serialize",
                          tolerance = 1e-2)
    expect_snapshot_value(res1$P2,
                          style = "serialize",
                          tolerance = 1e-2)

    res2 = contourmap.inla(data$result, data$stack, tag = "pred",
                           n.levels = 4, seed = data$seed,
                           max.threads = 1, ind = ind,
                           compute = list(F = FALSE, measures = c("P2","P1")),
                           method='QC')
    expect_snapshot_value(res2$P1,
                          style = "serialize",
                          tolerance = 1e-2)
    expect_snapshot_value(res2$P2,
                          style = "serialize",
                          tolerance = 1e-2)
    
    data2 <- testdata.inla.small(inla.mode = "compact")
    res3 = contourmap.inla(data2$result, data2$stack, tag = "pred",
                           n.levels = 4, seed = data2$seed,
                           max.threads = 1, ind = ind,
                           compute = list(F = FALSE, measures = c("P2","P1")),
                           method='EB')
    expect_equal(res1$P1, res3$P1, tolerance = 1e-2)
    expect_equal(res1$P2, res3$P2, tolerance = 1e-2)
    
    res4 = contourmap.inla(data2$result, data2$stack, tag = "pred",
                           n.levels = 4, seed = data2$seed,
                           max.threads = 1, ind = ind,
                           compute = list(F = FALSE, measures = c("P2","P1")),
                           method='QC')
    
    expect_equal(res2$P1, res4$P1, tolerance = 1e-2)
    expect_equal(res2$P2, res4$P2, tolerance = 6e-2)

})
