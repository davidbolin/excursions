## regions.inla.R
##
##   Copyright (C) 2026, David Bolin
##
##   This program is free software: you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation, either version 3 of the License, or
##   (at your option) any later version.
##
##   This program is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##
##   You should have received a copy of the GNU General Public License
##   along with this program.  If not, see <http://www.gnu.org/licenses/>.

#' Connected excursion regions for latent Gaussian models
#'
#' Connected excursion regions for latent Gaussian models fitted with `INLA`
#' or `inlabru`. See [excursions.regions()] for details on the regions.
#'
#' @param result.inla Result object from an `INLA` or `inlabru` call.
#' @param stack The stack object used in the INLA call.
#' @param name The name of the component for which to do the calculation. This
#' argument should only be used if a stack object is not provided, use the tag
#' argument otherwise. For `inlabru` results, use `name = "APredictor"`
#' together with `ind = bru_index(result, tag)`.
#' @param tag The tag of the component in the stack for which to do the
#' calculation. This argument should only be used if a stack object is
#' provided, use the name argument otherwise.
#' @param ind If only a part of a component should be used in the
#' calculations, this argument specifies the indices for that part.
#' @param method Method for handling the latent Gaussian structure:
#' \describe{
#' \item{'EB' }{Empirical Bayes}
#' \item{'QC' }{Quantile correction}
#' \item{'NI' }{Numerical integration}
#' \item{'NIQC' }{Numerical integration with quantile correction}
#' }
#' @param alpha Error probability for each region.
#' @param u Excursion level.
#' @param u.link If u.link is TRUE, `u` is assumed to be in the scale of the
#' data and is then transformed to the scale of the linear predictor (default
#' FALSE).
#' @param type Type of region, `'>'` for positive excursion regions and `'<'`
#' for negative excursion regions.
#' @param graph The neighbourhood graph of the nodes, given as a symmetric
#' sparse matrix where the non-zero off-diagonal elements are the edges, or
#' as an `fm_mesh_2d` object, in which case the vertex graph of the mesh is
#' used. The graph can either have one node for each of the nodes selected
#' by `ind`, in the same order, or one node for each node of the component.
#' The default is the graph of the non-zero elements of the precision
#' matrix, which for the linear predictor usually has no edges, so the graph
#' should then be provided.
#' @param n.iter Number or iterations in the MC sampler that is used for
#' approximating probabilities. The default value is 20000. If `size.tol` or `tol` is given, this is the maximal number of iterations.
#' @param max.regions The maximum number of regions to compute.
#' @param min.size The minimum number of nodes of a region.
#' @param growth How the regions are grown, see [excursions.regions()].
#' @param n.starts The maximum number of start points for growing a region in
#' each connected component, see [excursions.regions()].
#' @param min.prominence The minimum prominence of a local maximum for it to
#' be used as a start point, see [excursions.regions()].
#' @param verbose Set to TRUE for verbose mode (optional).
#' @param max.threads The number of threads that the program can use. The
#'   default, 0, uses the default number of threads of OpenMP.
#' @param compressed If INLA is run in compressed mode and a part of the
#' linear predictor is to be used, then only add the relevant part. Otherwise
#' the entire linear predictor is added internally (default TRUE).
#' @param seed Random seed (optional).
#' @param prune.ind If `TRUE` and `ind` is supplied, then the result object
#' is pruned to contain only the active nodes specified by `ind`, and the
#' regions are given as indices within `ind`.
#'
#' @param tol Target for the estimated errors of the joint probabilities of
#' the regions (optional), see [excursions.regions()].
#' @param size.tol Target for the estimated Monte Carlo error of the size of
#' each region, relative to the size, see [excursions.regions()].
#' @return `excursions.regions.inla` returns a list with the elements
#' \item{regions}{A list with the indices of the regions, largest first. The
#' indices refer to the nodes of the component, or to the positions in `ind`
#' if `prune.ind = TRUE`.}
#' \item{P}{The estimated joint excursion probability of each region.}
#' \item{P.err}{The Monte Carlo standard errors of `P`.}
#' \item{labels}{A vector with the region of each node, where 0 means that
#' the node is not in a region, and `NA` that it is not in `ind`.}
#' \item{E}{The largest region, as an indicator vector.}
#' \item{F}{The excursion functions of the regions, as a sparse matrix with
#' one column for each region, see [excursions.regions()].}
#' \item{rho}{Marginal excursion probabilities.}
#' \item{mean}{Posterior mean.}
#' \item{vars}{Marginal variances.}
#' \item{meta}{A list containing various information about the calculation.}
#' @export
#' @details The methods for handling the latent Gaussian structure are the
#' same as in [excursions.inla()], except that the `iNIQC` method is not
#' available. With the `NI` and `NIQC` methods, the joint probability of a
#' region is the mixture over the hyperparameter configurations of INLA, and
#' the regions are grown using the correlations of the configuration with
#' the largest posterior density.
#'
#' Models fitted with `inlabru` are handled in the same way as in
#' [excursions.inla()]. To compute regions for the linear predictor at a set
#' of locations, such as the nodes of a mesh, add a likelihood component with
#' `NA` observations at the locations and a tag, and use
#' `name = "APredictor"` and `ind = bru_index(result, tag)`. The graph is then
#' given for the locations, for example as the mesh.
#'
#' @note This function requires the `INLA` package, which is not a CRAN
#' package.  See <https://www.r-inla.org/download-install> for easy
#' installation instructions.
#' @author David Bolin \email{davidbolin@@gmail.com}
#' @seealso [excursions.regions()], [excursions.inla()]
#'
#' @examples
#' \dontrun{
#' if (require.nowarnings("INLA") && require.nowarnings("inlabru")) {
#'   ## Simulate data on a mesh
#'   x <- seq(from = 0, to = 10, length.out = 20)
#'   lattice <- fmesher::fm_lattice_2d(x = x, y = x)
#'   mesh <- fmesher::fm_rcdt_2d_inla(
#'     lattice = lattice, extend = FALSE, refine = FALSE
#'   )
#'   Q <- fmesher::fm_matern_precision(mesh, alpha = 2, rho = 3, sigma = 1)
#'   field <- fmesher::fm_sample(n = 1, Q = Q)
#'   obs.loc <- matrix(runif(200) * 10, 100, 2)
#'   y <- as.vector(fmesher::fm_basis(mesh, loc = obs.loc) %*% field) +
#'     rnorm(100) * 0.3
#'
#'   ## Fit the model with inlabru, with NA observations at the mesh nodes
#'   matern <- INLA::inla.spde2.pcmatern(mesh,
#'     prior.range = c(1, 0.5), prior.sigma = c(1, 0.5)
#'   )
#'   data <- data.frame(x1 = obs.loc[, 1], x2 = obs.loc[, 2], y = y)
#'   data.prd <- data.frame(x1 = mesh$loc[, 1], x2 = mesh$loc[, 2], y = NA)
#'   fit <- inlabru::bru(
#'     ~ Intercept(1) + field(cbind(x1, x2), model = matern),
#'     inlabru::bru_obs(y ~ ., family = "normal", data = data),
#'     inlabru::bru_obs(y ~ ., family = "normal", data = data.prd, tag = "prd"),
#'     options = list(control.compute = list(return.marginals.predictor = TRUE))
#'   )
#'
#'   ## Connected regions where the field exceeds 0
#'   res <- excursions.regions.inla(fit,
#'     name = "APredictor", ind = inlabru::bru_index(fit, "prd"),
#'     graph = mesh, alpha = 0.1, u = 0, type = ">", method = "QC",
#'     min.size = 5, prune.ind = TRUE
#'   )
#'   lengths(res$regions)
#' }
#' }
excursions.regions.inla <- function(result.inla,
                                    stack,
                                    name = NULL,
                                    tag = NULL,
                                    ind = NULL,
                                    method,
                                    alpha,
                                    u,
                                    u.link = FALSE,
                                    type,
                                    graph,
                                    n.iter = 20000,
                                    max.regions = Inf,
                                    min.size = 1,
                                    growth = c("bound", "rho"),
                                    n.starts = 10,
                                    min.prominence = 0,
                                    verbose = 0,
                                    max.threads = 0,
                                    compressed = TRUE,
                                    seed = NULL,
                                    prune.ind = FALSE,
                                    tol = NULL,
                                    size.tol = 0.001) {
  if (!requireNamespace("INLA", quietly = TRUE)) {
    stop("This function requires the INLA package (see www.r-inla.org/download-install)")
  }
  growth <- match.arg(growth)
  if (missing(result.inla)) {
    stop("Must supply INLA result object")
  }
  if (missing(method)) {
    cat("No method selected, using QC\n")
    method <- "QC"
  }
  if (!(method %in% c("EB", "QC", "NI", "NIQC"))) {
    stop("Method must be one of EB, QC, NI, NIQC")
  }
  if (missing(alpha)) {
    stop("Must specify error probability")
  }
  if (missing(u)) {
    stop("Must specify level u")
  }
  if (missing(type)) {
    stop("Must specify type of excursion")
  }
  if (!(type %in% c(">", "<"))) {
    stop("Only the types '>' and '<' are supported.")
  }
  if (!result.inla$.args$control.compute$config) {
    stop("INLA result must be calculated using control.compute$config=TRUE")
  }

  # Get indices for the component of interest in the configs
  tmp <- inla.output.indices(result.inla,
    name = name, stack = stack,
    tag = tag, compressed = compressed
  )
  ind.stack <- tmp$index
  if (tmp$result.updated) {
    result.inla <- tmp$result
    ind.stack.original <- tmp$index.original
  } else {
    ind.stack.original <- ind.stack
  }
  n <- length(result.inla$misc$configs$config[[1]]$mean)
  n.out <- length(ind.stack)
  ind.int <- seq_len(n.out)
  # ind is assumed to contain indices within the component of interest
  if (!is.null(ind)) {
    ind.int <- ind.int[ind]
    ind.stack <- ind.stack[ind]
    ind.stack.original <- ind.stack.original[ind]
  }
  ind <- ind.stack
  ind.original <- ind.stack.original

  # If u.link is TRUE, the limit is given in linear scale
  # then transform to the scale of the linear predictor
  u.t <- rho <- rep(0, n)
  if (u.link) {
    links <- result.inla$misc$linkfunctions$names[
      result.inla$misc$linkfunctions$link
    ]
    u.t[ind] <- sapply(ind, function(i) private.link.function(u, links[i]))
  } else {
    u.t <- u
  }

  if (verbose) {
    cat("Calculating marginal probabilities\n")
  }
  random.effect <- !is.null(name) && (name != "APredictor") &&
    (name != "Predictor")
  if (random.effect && is.null(result.inla$marginals.random)) {
    stop("INLA result must be calculated using return.marginals.random=TRUE if excursion sets to be calculated for a random effect of the model")
  }
  if (!random.effect && is.null(result.inla$marginals.linear.predictor)) {
    stop("INLA result must be calculated using return.marginals.predictor=TRUE if excursion sets are to be calculated for the linear predictor.")
  }
  if (random.effect) {
    rho.ind <- inla.get.marginal(ind.int,
      u = u, result = result.inla,
      effect.name = name, u.link = u.link, type = type
    )
  } else {
    rho.ind <- inla.get.marginal(ind.original,
      u = u, result = result.inla, u.link = u.link, type = type
    )
  }
  rho[ind] <- rho.ind

  ## The configurations, with weights given by the posterior densities
  mode.config <- private.get.config(result.inla, private.mode.config(result.inla))
  if (method == "EB" || method == "QC") {
    configs <- list(mode.config)
  } else {
    n.theta <- result.inla$misc$configs$nconfig
    configs <- lapply(seq_len(n.theta), function(i) {
      private.get.config(result.inla, i)
    })
    ## The configuration with the largest posterior density first, since the
    ## regions are grown using its correlations
    configs <- configs[order(-vapply(configs, function(x) x$lp, 0))]
  }
  weights <- exp(vapply(configs, function(x) x$lp, 0))
  configs <- lapply(configs, function(x) {
    list(mu = x$mu - u.t, Q = private.as.dgCMatrix(x$Q), vars = x$vars)
  })

  ## The graph over all nodes, with edges between the selected nodes
  if (missing(graph) || is.null(graph)) {
    graph <- configs[[1]]$Q[ind, ind, drop = FALSE]
  } else if (inherits(graph, c("fm_mesh_2d", "inla.mesh"))) {
    graph <- graph$graph$vv
  }
  if (all(dim(graph) == length(ind))) {
    graph.ind <- graph
  } else if (all(dim(graph) == n.out)) {
    graph.ind <- graph[ind.int, ind.int, drop = FALSE]
  } else {
    stop(paste0(
      "The graph must have one node for each of the selected nodes (",
      length(ind), ") or for each node of the component (", n.out, ")."
    ))
  }
  tr <- private.sparse.gettriplet(graph.ind)
  keep <- tr$x != 0 & tr$i != tr$j
  if (!any(keep) && length(ind) > 1) {
    warning("The graph has no edges between the selected nodes, so each region is a single node. Use the graph argument to give the neighbourhood structure, for example the mesh.")
  }
  G <- private.regions.graph(sparseMatrix(
    i = ind[tr$i[keep]], j = ind[tr$j[keep]], x = 1, dims = c(n, n)
  ), NULL, n)

  if (verbose) {
    cat("Calculating regions using the ", method, " method\n")
  }
  res <- private.regions(
    alpha = alpha, u = 0, configs = configs, weights = weights, type = type,
    qc = method %in% c("QC", "NIQC"), rho = rho, G = G, ind = ind,
    n.iter = n.iter, max.regions = max.regions, min.size = min.size,
    growth = growth, n.starts = n.starts, min.prominence = min.prominence,
    max.threads = max.threads, seed = seed, verbose = verbose, tol = tol,
    size.tol = size.tol
  )

  ## Indices in the component, or in ind if the result is pruned
  if (prune.ind) {
    pos <- seq_along(ind)
    n.res <- length(ind)
  } else {
    pos <- ind.int
    n.res <- n.out
  }
  regions <- lapply(res$regions, function(R) sort(pos[match(R, ind)]))
  labels <- E <- mu.out <- vars.out <- rho.out <- rep(NA, n.res)
  labels[pos] <- res$labels[ind]
  E[pos] <- res$E[ind]
  F.out <- Matrix(0, n.res, ncol(res$F), sparse = TRUE)
  F.out[pos, ] <- res$F[ind, , drop = FALSE]
  rho.out[pos] <- rho.ind
  mu.out[pos] <- mode.config$mu[ind]
  vars.out[pos] <- mode.config$vars[ind]

  output <- list(
    regions = regions,
    P = res$P,
    P.err = res$P.err,
    labels = labels,
    E = E,
    F = private.as.dgCMatrix(F.out),
    rho = rho.out,
    mean = mu.out,
    vars = vars.out,
    meta = list(
      calculation = "regions",
      type = type,
      level = u,
      level.link = u.link,
      alpha = alpha,
      n.iter = n.iter,
      tol = tol,
      size.tol = size.tol,
      method = method,
      growth = growth,
      n.starts = n.starts,
      min.prominence = min.prominence,
      min.size = min.size,
      max.regions = max.regions,
      ind = if (prune.ind) NULL else ind.int,
      call = match.call()
    )
  )
  class(output) <- "excurobj"
  output
}
