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
#' \code{contourmap} is used for calculating contour maps and quality measures for contour maps for Gaussian models.
#'
#' @param mu Expectation vector.
#' @param Q Precision matrix.
#' @param vars Precomputed marginal variances (optional).
#' @param n.levels Number of levels in contour map.
#' @param ind Indices of the nodes that should be analyzed (optional).
#' @param levels Levels to use in contour map.
#' @param type Type of contour map. One of:
#' \describe{
#'      \item{'standard' }{Equidistant levels between smallest and largest value of the posterior mean (default).}
#'      \item{'pretty' }{Equally spaced 'round' values which cover the range of the values in the posterior mean.}
#'      \item{'equalarea' }{Levels such that different spatial regions are approximately equal in size.}
#'      \item{'P0-optimal' }{Levels chosen to maximize the P0 measure.}
#'      \item{'P1-optimal' }{Levels chosen to maximize the P1 measure.}
#'      \item{'P2-optimal' }{Levels chosen to maximize the P2 measure.}
#' }
#' @param compute A list with quality indices to compute
#' \describe{
#'      \item{'F': }{TRUE/FALSE indicating whether the contour map function should be computed (default TRUE).}
#'      \item{'measures': }{A list with the quality measures to compute ("P0", "P1", "P2") or corresponding bounds based only on the marginal probabilities ("P0-bound", "P1-bound", "P2-bound").}
#'      }
#' @param use.marginals Only marginal distributions are used when finding P-optimal maps (default TRUE).
#' @param alpha Maximal error probability in contour map function (default=1).
#' @param F.limit The limit value for the computation of the F function. F is set to NA for all nodes where F<1-F.limit. Default is F.limit = \code{alpha}.
#' @param n.iter Number or iterations in the MC sampler that is used for calculating the quantities in \code{compute}. The default value is 10000.
#' @param verbose Set to TRUE for verbose mode (optional).
#' @param max.threads Decides the number of threads the program can use. Set to 0 for using the maximum number of threads allowed by the system (default).
#' @param seed Random seed (optional).
#'
#' @return \code{contourmap} returns an object of class "excurobj" with the following elements
#'     \item{u }{Contour levels used in the contour map.}
#'     \item{n.levels }{The number of contours used.}
#'     \item{u.e }{The values associated with the level sets G_k.}
#'     \item{G }{A vector which shows which of the level sets G_k each node belongs to.}
#'     \item{map }{Representation of the contour map with map[i]=u.e[k] if i is in G_k.}
#'     \item{F }{The contour map function (if computed).}
#'     \item{M }{Contour avoiding sets (if \code{F} is computed). \eqn{M=-1} for all non-significant nodes and  \eqn{M=k} for nodes that belong to \eqn{M_k}.}
#'     \item{P0/P1/P2 }{Calculated quality measures (if computed).}
#'     \item{P0bound/P1bound/P2bound }{Calculated upper bounds quality measures (if computed).}
#'     \item{meta }{A list containing various information about the calculation.}
#' @export
#' @details
#' The Gaussian model is specified using the mean \code{mu} and the precision matrix
#' \code{Q}. The contour map is then computed for the mean, using either the contour
#' levels specified in \code{levels}, or \code{n.levels} contours that are placed according
#' to the argument \code{type}.
#'
#' A number of quality measures can be computed based based on the specified contour map
#' and the Gaussian distribution. What should be computed is specified using the
#' \code{compute} argument. For details on these quanties, see the reference below.
#'
#' @author David Bolin \email{davidbolin@@gmail.com}
#' @references Bolin, D. and Lindgren, F. (2017) \emph{Quantifying the uncertainty of contour maps}, Journal of Computational and Graphical Statistics, vol 26, no 3, pp 513-524.
#'
#' Bolin, D. and Lindgren, F. (2018), \emph{Calculating Probabilistic Excursion Sets and Related Quantities Using excursions}, Journal of Statistical Software, vol 86, no 1, pp 1-20.
#' @seealso \code{\link{contourmap.inla}}, \code{\link{contourmap.mc}}, \code{\link{contourmap.colors}}
#' @examples
#' n <- 10
#' Q <- Matrix(toeplitz(c(1, -0.5, rep(0, n - 2))))
#' mu <- seq(-5, 5, length = n)
#' lp <- contourmap(mu, Q,
#'   n.levels = 2,
#'   compute = list(F = FALSE, measures = c("P1", "P2")),
#'   max.threads = 1
#' )
#' # Plot the contourmap
#' plot(lp$map)
#' # Display the quality measures
#' cat(c(lp$P1, lp$P2))
contourmap <- function(mu,
                       Q,
                       vars,
                       n.levels,
                       ind,
                       levels,
                       type = c(
                         "standard",
                         "pretty",
                         "equalarea",
                         "P0-optimal",
                         "P1-optimal",
                         "P2-optimal"
                       ),
                       compute = list(F = TRUE, measures = NULL),
                       use.marginals = TRUE,
                       alpha,
                       F.limit,
                       n.iter = 10000,
                       verbose = FALSE,
                       max.threads = 0,
                       seed = NULL) {
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
    if (missing(levels) || is.null(levels)) {
      stop("Must supply levels or n.levels")
    } else {
      n.levels <- length(levels)
    }
  }

  if (!missing(mu)) {
    mu <- private.as.vector(mu)
  }

  if (!missing(vars)) {
    vars <- private.as.vector(vars)
  }

  if (!missing(ind)) {
    ind <- private.as.vector(ind)
  }

  if (!missing(Q)) {
    Q <- private.as.dgCMatrix(Q)
  }


  measure <- NULL
  if (!is.null(compute$measures)) {
    measure <- match.arg(compute$measures,
      c("P0", "P1", "P2", "P0-bound", "P1-bound", "P2-bound"),
      several.ok = TRUE
    )
  }

  if (type == "standard") {
    if (verbose) cat("Creating contour map\n")
    lp <- excursions.levelplot(
      mu = mu, n.levels = n.levels, ind = ind,
      levels = levels, equal.area = FALSE
    )
  } else if (type == "pretty") {
    if (verbose) cat("Creating pretty contour map\n")
    lp <- excursions.levelplot(
      mu = mu, n.levels = n.levels, ind = ind,
      levels = levels, equal.area = FALSE, pretty.cm = TRUE
    )
    n.levels <- lp$n.levels
  } else if (type == "equalarea") {
    if (verbose) cat("Creating equal area contour map\n")
    lp <- excursions.levelplot(
      mu = mu, n.levels = n.levels, ind = ind,
      levels = levels, equal.area = TRUE
    )
  } else if (type == "P0-optimal" || type == "P1-optimal" || type == "P2-optimal") {
    if (!missing(levels)) {
      warning("Not using supplied levels for optimal contour map\n")
      if (!missing(n.levels)) {
        if (n.levels != length(levels)) {
          warning("n.levels != length(levels), using n.levels\n")
        }
      } else {
        n.levels <- length(levels)
      }
    }
    if (missing(vars) && missing(Q)) {
      stop("Variances must be supplied when creating optimal contour map")
    } else if (missing(vars)) {
      vars <- excursions.variances(Q = Q, max.threads = max.threads)
    }
    if (use.marginals == TRUE) {
      if (missing(Q)) {
        stop("The precision matrix must be supplied unless marginals are used")
      }
    }

    if (type == "P0-optimal") {
      if (verbose) cat("Creating P0-optimal contour map\n")
      opt.measure <- 0
    } else if (type == "P1-optimal") {
      if (verbose) cat("Creating P1-optimal contour map\n")
      opt.measure <- 1
    } else if (type == "P2-optimal") {
      if (verbose) cat("Creating P2-optimal contour map\n")
      opt.measure <- 2
    }

    lp <- excursions.opt.levelplot(
      mu = mu, vars = vars, Q = Q,
      n.levels = n.levels, measure = opt.measure,
      use.marginals = use.marginals, ind = ind
    )
  }

  F.calculated <- FALSE
  if (!is.null(measure)) {
    if (missing(Q)) {
      stop("precision matrix must be supplied if measure should be calculated")
    }

    for (i in seq_along(measure)) {
      if (measure[i] == "P1") {
        if (n.levels > 1) {
          if (verbose) cat("Calculating P1-measure\n")
          tmp <- Pmeasure(lp = lp, mu = mu, Q = Q, ind = ind, type = 1, 
                          seed = seed, n.iter = n.iter, max.threads = max.threads)
          lp$P1 <- tmp$P
          lp$P1.error <- tmp$E
        } else {
          lp$P1 <- 1
          lp$P1.error <- 0
        }
      } else if (measure[i] == "P2") {
        if (verbose) cat("Calculating P2-measure\n")
        tmp <- Pmeasure(lp = lp, mu = mu, Q = Q, ind = ind, type = 2, 
                        seed = seed, n.iter = n.iter, max.threads = max.threads)
        lp$P2 <- tmp$P
        lp$P2.error <- tmp$E
      } else if (measure[i] == "P0") {
        if (verbose) cat("Calculating P0-measure and contour map function\n")

        p <- contourfunction(
          lp = lp, mu = mu, Q = Q, vars = vars, ind = ind,
          alpha = alpha, F.limit = F.limit,
          n.iter = n.iter, max.threads = max.threads,
          seed = seed, verbose = verbose
        )
        F.calculated <- TRUE
      } else if (measure[i] == "P0-bound") {
        if (missing(vars)) {
          vars <- excursions.variances(Q = Q, max.threads = max.threads)
        }
        lp$P0.bound <- Pmeasure.bound(lp = lp, mu = mu, vars, type = 0, ind = ind)
      } else if (measure[i] == "P1-bound") {
        if (missing(vars)) {
          vars <- excursions.variances(Q = Q, max.threads = max.threads)
        }
        lp$P1.bound <- Pmeasure.bound(lp = lp, mu = mu, vars, type = 1, ind = ind)
      } else if (measure[i] == "P2-bound") {
        if (missing(vars)) {
          vars <- excursions.variances(Q = Q, max.threads = max.threads)
        }
        lp$P2.bound <- Pmeasure.bound(lp = lp, mu = mu, vars, type = 2, ind = ind)
      }
    }
  }
  if (!F.calculated) {
    if (is.null(compute$F) || compute$F) {
      if (verbose) cat("Calculating contour map function\n")
      p <- contourfunction(
        lp = lp, mu = mu, Q = Q, vars = vars, ind = ind,
        alpha = alpha, F.limit = F.limit,
        n.iter = n.iter, max.threads = max.threads,
        seed = seed, verbose = verbose
      )
      F.calculated <- TRUE
    }
  }

  if (missing(ind) || is.null(ind)) {
    ind <- seq_len(length(mu))
  } else if (is.logical(ind)) {
    ind <- which(ind)
  }

  if (F.calculated) {
    lp$P0 <- mean(p$F[ind])
    lp$F <- p$F
    lp$E <- p$E
    lp$M <- p$M
    lp$rho <- p$rho
    # } else {
    # lp$E <- NULL
  }
  lp$meta <- list(
    calculation = "contourmap",
    F.limit = F.limit,
    F.computed = compute$F,
    alpha = alpha,
    levels = lp$u,
    type = "!=",
    contourmap.type = type,
    n.iter = n.iter,
    mu.range = range(mu[ind]),
    mu = mu[ind],
    ind = ind,
    call = match.call()
  )
  class(lp) <- "excurobj"
  return(lp)
}
