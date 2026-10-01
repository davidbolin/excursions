test_that("excursions.regions.inla, methods", {
  skip_on_cran()
  local_exc_safe_inla()

  data <- testdata.inla()
  graph <- bandSparse(data$n, k = 1, symmetric = TRUE)
  for (method in c("EB", "QC", "NI", "NIQC")) {
    res <- excursions.regions.inla(data$result, data$stack,
      tag = "pred", method = method, alpha = 0.1, u = 0, type = ">",
      graph = graph, seed = data$seed, max.threads = 1
    )
    expect_true(length(res$regions) >= 1)
    for (R in res$regions) {
      expect_true(all(diff(R) == 1))
    }
    expect_false(anyDuplicated(unlist(res$regions)) > 0)
    expect_true(all(res$P >= 0.9))
    expect_equal(which(res$labels == 1), res$regions[[1]])
    expect_true(all(res$rho[unlist(res$regions)] >= 0.9))
  }
})

test_that("excursions.regions.inla, ind and graph", {
  skip_on_cran()
  local_exc_safe_inla()

  data <- testdata.inla()
  ind <- 3:9
  graph <- bandSparse(data$n, k = 1, symmetric = TRUE)
  res <- excursions.regions.inla(data$result, data$stack,
    tag = "pred", ind = ind, method = "QC", alpha = 0.1, u = 0, type = ">",
    graph = graph, seed = data$seed, max.threads = 1
  )
  expect_true(all(unlist(res$regions) %in% ind))
  expect_true(all(is.na(res$labels[-ind])))

  ## A graph for the selected nodes gives the same result as a graph for
  ## the component, and pruning gives the positions within ind
  res.p <- excursions.regions.inla(data$result, data$stack,
    tag = "pred", ind = ind, method = "QC", alpha = 0.1, u = 0, type = ">",
    graph = graph[ind, ind], seed = data$seed, max.threads = 1,
    prune.ind = TRUE
  )
  expect_equal(res.p$regions, lapply(res$regions, function(R) match(R, ind)))
  expect_equal(res.p$P, res$P)
  expect_equal(res.p$labels, res$labels[ind])
  expect_equal(dim(res$F), c(data$n, length(res$regions)))
  expect_equal(as.matrix(res.p$F), as.matrix(res$F[ind, , drop = FALSE]))
  expect_s3_class(res, "excurobj")

  expect_error(excursions.regions.inla(data$result, data$stack,
    tag = "pred", ind = ind, method = "QC", alpha = 0.1, u = 0, type = ">",
    graph = graph[1:5, 1:5]
  ))
  expect_error(excursions.regions.inla(data$result, data$stack,
    tag = "pred", method = "iNIQC", alpha = 0.1, u = 0, type = ">",
    graph = graph
  ))
  expect_warning(excursions.regions.inla(data$result, data$stack,
    tag = "pred", method = "QC", alpha = 0.1, u = 0, type = ">",
    max.regions = 1, seed = data$seed, max.threads = 1
  ))
})

test_that("excursions.regions.inla, inlabru", {
  skip_on_cran()
  skip_if_not_installed("inlabru")
  local_exc_safe_inla()

  set.seed(4)
  x <- seq(0, 10, length.out = 10)
  mesh <- fmesher::fm_rcdt_2d_inla(
    lattice = fmesher::fm_lattice_2d(x = x, y = x),
    extend = FALSE, refine = FALSE
  )
  Q <- fmesher::fm_matern_precision(mesh, alpha = 2, rho = 4, sigma = 1)
  field <- fmesher::fm_sample(n = 1, Q = Q)
  obs <- matrix(runif(160) * 10, 80, 2)
  y <- 0.5 + as.vector(fmesher::fm_basis(mesh, loc = obs) %*% field) +
    rnorm(80) * 0.3
  matern <- INLA::inla.spde2.pcmatern(mesh,
    prior.range = c(1, 0.5), prior.sigma = c(1, 0.5)
  )
  data <- data.frame(x1 = obs[, 1], x2 = obs[, 2], y = y)
  prd <- data.frame(x1 = mesh$loc[, 1], x2 = mesh$loc[, 2], y = NA)
  fit <- inlabru::bru(~ Intercept(1) + field(cbind(x1, x2), model = matern),
    inlabru::bru_obs(y ~ ., family = "normal", data = data),
    inlabru::bru_obs(y ~ ., family = "normal", data = prd, tag = "prd"),
    options = list(
      control.compute = list(return.marginals.predictor = TRUE),
      num.threads = "1:1"
    )
  )
  res <- excursions.regions.inla(fit,
    name = "APredictor", ind = inlabru::bru_index(fit, "prd"),
    graph = mesh, alpha = 0.1, u = 0, type = ">", method = "QC",
    prune.ind = TRUE, seed = 1, max.threads = 1
  )
  expect_length(res$labels, mesh$n)
  expect_true(all(res$P >= 0.9))
  G <- private.regions.graph(mesh, NULL, mesh$n)
  for (R in res$regions) {
    mask <- integer(mesh$n)
    mask[R] <- 1L
    lab <- .Call("regions_components", G@p, G@i, mask, PACKAGE = "excursions")
    expect_equal(max(lab), 1L)
  }
})
