test_that("excursions.inla, test ind", {
  skip_on_cran()
  local_exc_safe_inla()

  data <- testdata.inla()
  ind1 <- c(1, 2, 3, 4)
  ind2 <- c(4, 3, 2, 1)
  ind3 <- rep(FALSE, data$n)
  ind3[1:4] <- TRUE

  res1 <- excursions.inla(data$result, data$stack,
    ind = ind1, method = "QC",
    tag = "pred", u = 0, type = ">", seed = data$seed,
    max.threads = 1
  )
  res2 <- excursions.inla(data$result, data$stack,
    ind = ind2, method = "QC",
    tag = "pred", u = 0, type = ">", seed = data$seed,
    max.threads = 1
  )
  res3 <- excursions.inla(data$result, data$stack,
    ind = ind3, method = "QC",
    tag = "pred", u = 0, type = ">", seed = data$seed,
    max.threads = 1
  )

  expect_equal(res1$F, res2$F, tolerance = 1e-7)
  expect_equal(res2$F, res3$F, tolerance = 1e-7)
})

test_that("excursions.inla, compact mode", {
  skip_on_cran()
  local_exc_safe_inla()

  data1 <- testdata.inla(inla.mode = "classic")
  data2 <- testdata.inla(inla.mode = "compact")

  res1 <- excursions.inla(data1$result,
    name = "ar", method = "QC",
    u = 0, type = ">", seed = data1$seed,
    max.threads = 1
  )
  res2 <- excursions.inla(data2$result,
    name = "ar", method = "QC",
    u = 0, type = ">", seed = data2$seed, compressed = FALSE,
    max.threads = 1
  )
  res3 <- excursions.inla(data2$result,
    name = "ar", method = "QC",
    u = 0, type = ">", seed = data2$seed,
    max.threads = 1
  )

  # The inla estimates for different inla.mode will be different,
  # but should be similar
  expect_equal(res1$F, res2$F, tolerance = 1e-2)
  expect_equal(res2$F, res3$F, tolerance = 1e-4)
})

test_that("excursions.inla, compact mode, indexing", {
  skip_on_cran()
  local_exc_safe_inla()

  data1 <- testdata.inla(inla.mode = "classic")
  data2 <- testdata.inla(inla.mode = "compact")

  ind1 <- c(1, 2, 3, 4)
  ind2 <- c(4, 3, 2, 1)
  ind3 <- rep(FALSE, data1$n)
  ind3[1:4] <- TRUE

  res1 <- excursions.inla(data1$result,
    stack = data1$stack, tag = "pred",
    method = "QC", u = 0, type = ">", seed = data1$seed,
    ind = ind1,
    max.threads = 1
  )
  res2 <- excursions.inla(data1$result,
    stack = data1$stack, tag = "pred",
    method = "QC", u = 0, type = ">", seed = data1$seed,
    ind = ind2,
    max.threads = 1
  )
  res3 <- excursions.inla(data1$result,
    stack = data1$stack, tag = "pred",
    method = "QC", u = 0, type = ">", seed = data1$seed,
    ind = ind3,
    max.threads = 1
  )
  expect_equal(res1$F, res2$F, tolerance = 1e-2)
  expect_equal(res2$F, res3$F, tolerance = 1e-4)

  res4 <- excursions.inla(data2$result,
    stack = data2$stack, tag = "pred",
    method = "QC", u = 0, type = ">", seed = data1$seed,
    ind = ind1, compressed = FALSE,
    max.threads = 1
  )
  res5 <- excursions.inla(data2$result,
    stack = data2$stack, tag = "pred",
    method = "QC", u = 0, type = ">", seed = data1$seed,
    ind = ind2, compressed = FALSE,
    max.threads = 1
  )
  res6 <- excursions.inla(data2$result,
    stack = data2$stack, tag = "pred",
    method = "QC", u = 0, type = ">", seed = data1$seed,
    ind = ind3, compressed = FALSE,
    max.threads = 1
  )

  expect_equal(res3$F[1:4], res4$F[1:4], tolerance = 5e-2)
  expect_equal(res4$F, res5$F, tolerance = 1e-4)
  expect_equal(res5$F, res6$F, tolerance = 1e-4)

  res7 <- excursions.inla(data2$result,
    stack = data2$stack, tag = "pred",
    method = "QC", u = 0, type = ">", seed = data1$seed,
    ind = ind1, compressed = TRUE,
    max.threads = 1
  )
  res8 <- excursions.inla(data2$result,
    stack = data2$stack, tag = "pred",
    method = "QC", u = 0, type = ">", seed = data1$seed,
    ind = ind2, compressed = TRUE,
    max.threads = 1
  )
  res9 <- excursions.inla(data2$result,
    stack = data2$stack, tag = "pred",
    method = "QC", u = 0, type = ">", seed = data1$seed,
    ind = ind3, compressed = TRUE,
    max.threads = 1
  )

  expect_equal(res6$F[1:4], res7$F[1:4], tolerance = 5e-2)
  expect_equal(res7$F, res8$F, tolerance = 1e-4)
  expect_equal(res8$F, res9$F, tolerance = 1e-4)
})
