## excursions.inla.R
##
##   Copyright (C) 2012, 2013, 2014 David Bolin, Finn Lindgren
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

#' Contour maps and contour map quality measures for latent Gaussian models
#'
#' An interface to the \code{contourmap} function for latent Gaussian models
#' calculated using the INLA method.
#'
#' @param result.inla Result object from INLA call.
#' @param stack The stack object used in the INLA call.
#' @param name The name of the component for which to do the calculation. This
#' argument should only be used if a stack object is not provided, use the tag
#' argument otherwise.
#' @param tag The tag of the component in the stack for which to do the
#' calculation. This argument should only be used if a stack object is provided,
#' use the name argument otherwise.
#' @param method Method for handeling the latent Gaussian structure. Currently
#' only Empirical Bayes (EB) and Quantile corrections (QC) are supported.
#' @param n.levels Number of levels in contour map.
#' @param type Type of contour map. One of:
#'  \describe{
#'      \item{'standard' }{Equidistant levels between smallest and largest value
#'      of the posterior mean (default).}
#'      \item{'pretty' }{Equally spaced 'round' values which cover the range of
#'      the values in the posterior mean.}
#'      \item{'equalarea' }{Levels such that different spatial regions are
#'      approximately equal in size.}
#'      }
#' @param compute A list with quality indices to compute
#' \describe{
#'      \item{'F': }{TRUE/FALSE indicating whether the contour map function
#'      should be computed (default TRUE)}
#'      \item{'measures': }{A list with the quality measures to compute ("P0",
#'      "P1", "P2") or corresponding bounds based only on the marginal
#'      probabilities ("P0-bound", "P1-bound", "P2-bound")}
#'      }
#' @param alpha Maximal error probability in contour map function (default=1)
#' @param F.limit The limit value for the computation of the F function. F is
#' set to NA for all nodes where F<1-F.limit. Default is F.limit = \code{alpha}.
#' @param n.iter Number or iterations in the MC sampler that is used for
#' calculating the quantities in \code{compute}. The default value is 10000.
#' @param verbose Set to TRUE for verbose mode (optional)
#' @param max.threads Decides the number of threads the program can use. Set to
#' 0 for using the maximum number of threads allowed by the system (default).
#' @param seed Random seed (optional).
#' @param compressed If INLA is run in compressed mode and a part of the linear
#' predictor is to be used, then only add the relevant part. Otherwise the
#' entire linear predictor is added internally (default TRUE).
#' @param ind If only a part of a component should be used in the calculations,
#' this argument specifies the indices for that part (optional).
#' @param ... Additional arguments to the contour map function. See the
#' documentation for \code{contourmap} for details.
#'
#' @return \code{contourmap.inla} returns an object of class "excurobj" with the
#' same elements as returned by \code{contourmap}.
#' @note This function requires the \code{INLA} package, which is not a CRAN
#' package.  See \url{https://www.r-inla.org/download-install} for easy
#' installation instructions.
#' @author David Bolin \email{davidbolin@@gmail.com}
#' @details
#' The INLA approximation of the quantity of interest is in general a weighted
#' sum of Gaussian distributions with different parameters. If
#' \code{method = 'EB'} is used, then the contour map is computed for the mean
#' of the component in the weighted sum that has parameters with the highest
#' likelihood. If on the other hand \code{method='QC'}, then the contour map is
#' computed for the posterior mean reported by INLA. If the EB method also is
#' used in INLA, then this reported posterior mean is equal to the mean of the
#' component with the highest likelihood. Therefore, \code{method='EB'} is
#' appropriate if the EB method also is used in INLA, but \code{method='QC'}
#' should be used in general.
#'
#' The \code{n.levels} contours in the contour map are are placed according
#' to the argument \code{type}. A number of quality measures can be computed
#' based based on the specified contour map and the distribution of the
#' component of interest. What should be computed is specified using the
#' \code{compute} argument. For details on these quanties, see the reference
#' below.
#' @references Bolin, D. and Lindgren, F. (2017) \emph{Quantifying the
#' uncertainty of contour maps}, Journal of Computational and Graphical
#' Statistics, 26:3, 513-524.
#'
#' Bolin, D. and Lindgren, F. (2018), \emph{Calculating Probabilistic Excursion
#' Sets and Related Quantities Using excursions}, Journal of Statistical
#' Software, vol 86, no 1, pp 1-20.
#' @export
#' @seealso \code{\link{contourmap}}, \code{\link{contourmap.mc}},
#' \code{\link{contourmap.colors}}
#' @examples
#' \dontrun{
#' if (require.nowarnings("INLA")) {
#'   # Generate mesh and SPDE model
#'   n.lattice <- 10 # increase for more interesting, but slower, examples
#'   x <- seq(from = 0, to = 10, length.out = n.lattice)
#'   lattice <- fmesher::fm_lattice_2d(x = x, y = x)
#'   mesh <- fmesher::fm_rcdt_2d_inla(lattice = lattice, extend = FALSE, refine = FALSE)
#'   spde <- inla.spde2.matern(mesh, alpha = 2)
#'
#'   # Generate an artificial sample
#'   sigma2.e <- 0.01
#'   n.obs <- 100
#'   obs.loc <- cbind(
#'     runif(n.obs) * diff(range(x)) + min(x),
#'     runif(n.obs) * diff(range(x)) + min(x)
#'   )
#'   Q <- inla.spde2.precision(spde, theta = c(log(sqrt(0.5)), log(sqrt(1))))
#'   x <- inla.qsample(Q = Q)
#'   A <- fmesher::fm_basis(mesh = mesh, loc = obs.loc)
#'   Y <- as.vector(A %*% x + rnorm(n.obs) * sqrt(sigma2.e))
#'
#'   ## Estimate the parameters using INLA
#'   mesh.index <- inla.spde.make.index(name = "field", n.spde = spde$n.spde)
#'   ef <- list(c(mesh.index, list(Intercept = 1)))
#'
#'   s.obs <- inla.stack(data = list(y = Y), A = list(A), effects = ef, tag = "obs")
#'   s.pre <- inla.stack(data = list(y = NA), A = list(1), effects = ef, tag = "pred")
#'   stack <- inla.stack(s.obs, s.pre)
#'   formula <- y ~ -1 + Intercept + f(field, model = spde)
#'   result <- inla(
#'     formula = formula, family = "normal", data = inla.stack.data(stack),
#'     control.predictor = list(
#'       A = inla.stack.A(stack),
#'       compute = TRUE
#'     ),
#'     control.compute = list(
#'       config = TRUE,
#'       return.marginals.predictor = TRUE
#'     ),
#'     num.threads = 1
#'   )
#'
#'   ## Calculate contour map with two levels
#'   map <- contourmap.inla(result,
#'     stack = stack, tag = "pred",
#'     n.levels = 2, alpha = 0.1, F.limit = 0.1,
#'     max.threads = 1
#'   )
#'
#'   ## Plot the results
#'   cols <- contourmap.colors(map,
#'     col = heat.colors(100, 1),
#'     credible.col = grey(0.5, 1)
#'   )
#'   image(matrix(map$M[mesh$idx$lattice], n.lattice, n.lattice), col = cols)
#' }
#' }
contourmap.inla <- function(result.inla,
                            stack,
                            name = NULL,
                            tag = NULL,
                            method = "QC",
                            n.levels,
                            type = c(
                              "standard",
                              "pretty",
                              "equalarea"
                            ),
                            compute = list(F = TRUE, measures = NULL),
                            alpha,
                            F.limit,
                            n.iter = 10000,
                            verbose = FALSE,
                            max.threads = 0,
                            compressed = TRUE,
                            seed = NULL,
                            ind, ...) {
  if (!requireNamespace("INLA", quietly = TRUE)) {
    stop("This function requires the INLA package (see www.r-inla.org/download-install)")
  }
  if (missing(result.inla)) {
    stop("Must supply INLA result object")
  }

  type <- match.arg(type)

  if (missing(alpha) || is.null(alpha)) {
    alpha <- 0.1
  }
  if (missing(F.limit)) {
    F.limit <- 0.99
  } else {
    F.limit <- max(alpha, F.limit)
  }
  if (missing(n.levels) || is.null(n.levels)) {
    stop("Must supply n.levels")
  }
  measure <- NULL
  if (!is.null(compute$measures)) {
    measure <- match.arg(compute$measures,
      c("P0", "P1", "P2", "P0-bound", "P1-bound", "P2-bound"),
      several.ok = TRUE
    )
  }

  if (compute$F) { # compute P0 measure if F is computed anyway
    if (is.null(measure)) {
      measure <- c("P0")
    } else if (("P0" %in% measure) == FALSE) {
      measure <- c(measure, "P0")
    }
  }
  if (method != "EB" && method != "QC") {
    stop("Currently only EB and QC methods are implemented")
  }

  qc <- FALSE
  if (method == "QC") {
    qc <- TRUE
  }
  # compute indices, here ind will contain the indices that are used to extract the
  # relevant part from the configs, ind.int is the index vector for extracting marginal
  # distributions for random effects, and indices is a logical version of ind
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

  n.out <- length(ind.stack)
  ind.int <- seq_len(n.out)
  if (!missing(ind) && !is.null(ind)) {
    ind.int <- ind.int[ind]
    ind.stack <- ind.stack[ind]
    ind.stack.original <- ind.stack.original[ind]
  }
  ind <- ind.stack
  ind.original <- ind.stack.original

  for (i in 1:result.inla$misc$configs$nconfig) {
    config <- private.get.config(result.inla, i)
    if (config$lp == 0) {
      break
    }
  }
  indices <- rep(FALSE, length(config$mu))
  indices[ind] <- TRUE

  # Are we interested in a random effect?
  random.effect <- FALSE
  if (!missing(name) && !is.null(name) && name != "APredictor" && name != "Predictor") {
    random.effect <- TRUE
    if (is.null(result.inla$marginals.random)) {
      stop("INLA result must be calculated using return.marginals.random=TRUE if P measures are to be calculated for a random effect of the model")
    }
  }

  if (!random.effect && is.null(result.inla$marginals.linear.predictor)) {
    stop("INLA result must be calculated using return.marginals.predictor=TRUE if P measures are to be calculated for the linear predictor.")
  }

  # If method=EB, we use the mean of the configuration at the mode for the contourmap
  # If method=QC, we instead use the total mean.
  if (method == "EB") {
    mu <- config$mu
  } else {
    mu <- matrix(0, length(config$mu), 1)
    if (random.effect) {
      mu[ind] <- result.inla$summary.random[[name]][ind.int]
    } else {
      mu[ind] <- result.inla$summary.linear.predictor$mean[ind]
    }
  }

  # Compute the contour map
  cm <- contourmap(
    mu = config$mu,
    Q = config$Q,
    ind = ind,
    compute = list(F = FALSE, measures = NULL),
    n.levels = n.levels, ...
  )

  # compute measures
  F.computed <- FALSE
  if (!is.null(measure)) {
    rho <- matrix(0, length(config$mu), 2)
    for (i in seq_along(measure)) {
      # compute limits
      limits <- excursions.limits(lp = cm, mu = config$mu, measure = measure[i])
      # if QC, update limits
      if (method == "QC") {
        if (random.effect) {
          rho.ind <- sapply(seq_along(ind), function(j) {
            inla.get.marginal.int(ind.int[j],
              a = limits$a[ind[j]],
              b = limits$b[ind[j]],
              result = result.inla,
              effect.name = name
            )
          })
        } else {
          rho.ind <- sapply(seq_along(ind), function(j) {
            inla.get.marginal.int(ind.original[j],
              a = limits$a[ind[j]],
              b = limits$b[ind[j]],
              result = result.inla
            )
          })
        }
        rho[ind, ] <- t(rho.ind)
        limits$a <- config$mu + sqrt(config$vars) * qnorm(pmin(pmax(rho[, 1], 0), 1))
        limits$b <- config$mu + sqrt(config$vars) * qnorm(pmin(pmax(rho[, 2], 0), 1))
      } else {
        rho <- NULL
      }

      if (measure[i] == "P1") {
        if (n.levels > 1) {
          if (verbose) cat("Calculating P1-measure\n")
          tmp <- gaussint(
            mu = config$mu, Q = config$Q, a = limits$a,
            b = limits$b, ind = indices, use.reordering = "limits",
            n.iter = n.iter, seed = seed
          )
          cm$P1 <- tmp$P[1]
          cm$P1.error <- tmp$E[1]
        } else {
          cm$P1 <- 1
          cm$P1.error <- 0
        }
      } else if (measure[i] == "P2") {
        if (verbose) cat("Calculating P2-measure\n")
        tmp <- gaussint(
          mu = config$mu, Q = config$Q, a = limits$a,
          b = limits$b, ind = indices, use.reordering = "limits",
          max.threads = max.threads,
          n.iter = n.iter, seed = seed
        )
        cm$P2 <- tmp$P
        cm$P2.error <- tmp$E
      } else if (measure[i] == "P0") {
        if (verbose) cat("Calculating P0-measure and contour map function\n")
        p <- contourfunction(
          lp = cm, mu = config$mu, Q = config$Q,
          vars = config$vars, ind = ind, alpha = alpha,
          F.limit = F.limit, rho = rho, n.iter = n.iter,
          max.threads = max.threads, seed = seed,
          verbose = verbose, qc = qc
        )
        cm$P0 <- mean(p$F[ind])
        cm$F <- p$F
        cm$E <- p$E
        cm$M <- p$M
        cm$rho <- p$rho
        cm$meta$F.limit <- F.limit
        cm$meta$alpha <- alpha
        F.computed <- TRUE
      } else { # bounds
        if (method != "QC") {
          rho <- matrix(0, length(config$mu), 2)
          rho[ind, 1] <- pnorm(limits$a[ind], config$mu[ind], sqrt(config$vars[ind]))
          rho[ind, 2] <- pnorm(limits$b[ind], config$mu[ind], sqrt(config$vars[ind]))
        }
        if (measure[i] == "P0-bound") {
          cm$P0.bound <- mean(rho[ind, 2] - rho[ind, 1])
        } else if (measure[i] == "P1-bound") {
          cm$P1.bound <- min(rho[ind, 2] - rho[ind, 1])
        } else if (measure[i] == "P2-bound") {
          cm$P2.bound <- min(rho[ind, 2] - rho[ind, 1])
        }
      }
    }
  }
  # Fix output
  set.out <- rep(NA, n.out)
  if (F.computed) {
    set.out[ind.int] <- cm$E[ind]
    cm$E <- set.out
  }

  if (!is.null(cm$G)) {
    set.out[ind.int] <- cm$G[ind]
    cm$G <- set.out
  }

  if (!is.null(cm$F)) {
    F0 <- cm$F[ind]
    F0[is.na(F0)] <- 0
    cm$P0 <- mean(F0)
    set.out[ind.int] <- cm$F[ind]
    cm$F <- set.out
  }
  cm$meta$F.computed <- F.computed

  if (!is.null(cm$M)) {
    set.out[ind.int] <- cm$M[ind]
    cm$M <- set.out
  }

  if (!is.null(cm$rho)) {
    set.out[ind.int] <- cm$rho[ind]
    cm$rho <- set.out
  }


  if (!is.null(cm$map)) {
    set.out[ind.int] <- cm$map[ind]
    cm$map <- set.out
  }


  cm$meta$ind <- ind.int
  cm$meta$call <- match.call()
  return(cm)
}
