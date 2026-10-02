## regions.R
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

#' Connected Excursion Regions for Gaussian Random Fields
#'
#' `excursions.regions` computes connected (contiguous) regions where a
#' Gaussian random field jointly exceeds a level with high probability. The
#' first region is the largest connected region \eqn{D} that is found with
#' \eqn{P(X(s) > u, s \in D) \geq 1 - \alpha}{P(X(s) > u for all s in D) >= 1 - alpha}.
#' The nodes of the region are then removed, and the search is repeated for
#' the second largest region, and so on. Each region individually satisfies
#' the probability requirement, and the regions do not overlap.
#'
#' @param alpha Error probability for each region.
#' @param u Excursion level.
#' @param mu Expectation vector.
#' @param Q Precision matrix.
#' @param type Type of region, `'>'` for positive excursion regions and `'<'`
#' for negative excursion regions.
#' @param graph The neighbourhood graph that defines which nodes are
#' connected, given as a symmetric sparse matrix where the non-zero
#' off-diagonal elements are the edges, or as an `fm_mesh_2d` object, in
#' which case the vertex graph of the mesh is used. The default is the
#' graph of the non-zero elements of `Q`.
#' @param n.iter Number or iterations in the MC sampler that is used for
#' approximating probabilities. The default value is 20000. If `size.tol` or `tol` is given, this is the maximal number of iterations.
#' @param vars Precomputed marginal variances (optional).
#' @param rho Marginal excursion probabilities (optional).
#' @param method Method for handling the latent Gaussian structure:
#'  \describe{
#'       \item{'EB'}{Empirical Bayes (default)}
#'       \item{'QC'}{Quantile correction, rho must be provided if QC is used.}}
#' @param ind Indices of the nodes that can be included in the regions
#' (optional).
#' @param max.regions The maximum number of regions to compute. The default
#' computes all regions.
#' @param min.size The minimum number of nodes of a region. Regions with fewer
#' nodes are not computed.
#' @param growth How the regions are grown from the node with the largest
#' marginal probability:
#'  \describe{
#'       \item{'bound'}{The region grows by the neighbour that decreases a
#'       lower bound of the joint probability the least (default). This
#'       prefers nodes that are strongly correlated with the region.}
#'       \item{'rho'}{The region grows by the neighbour with the largest
#'       marginal excursion probability.}}
#' @param n.starts The maximum number of start points for growing a region in
#' each connected component. The start points are the local maxima of the
#' marginal excursion probabilities in the component, largest first, and the
#' largest region over the start points is kept. Local maxima in the largest
#' region found so far are skipped, and do not count towards `n.starts`.
#' @param min.prominence The minimum prominence of a local maximum for it to
#' be used as a start point, on the standardised scale
#' \eqn{(\mu - u)/\sigma}{(mu - u)/sigma}. The prominence is how far the
#' standardised margin must decrease from the maximum before a higher
#' maximum can be reached, so small values remove maxima that are caused by
#' noise in the mean. The largest maximum of each component is always used.
#' @param max.threads The number of threads that the program can use. The
#'   default, 0, uses the default number of threads of OpenMP.
#' @param seed Random seed (optional).
#' @param verbose Set to TRUE for verbose mode (optional).
#'
#' @param tol Target for the estimated errors of the joint probabilities of
#' the regions where they pass `1 - alpha` (optional). If `tol` is given, the
#' number of iterations is chosen adaptively, using at most `n.iter`
#' iterations, see [gaussint()]. It takes precedence over `size.tol`.
#' @param size.tol Target for the estimated Monte Carlo error of the size of
#' each region, relative to the size, when it is grown, see
#' [excursions()]. The default is 0.001. It is not used if `tol` is given.
#' @return `excursions.regions` returns a list with the elements
#' \item{regions}{A list with the node indices of the regions, largest first.}
#' \item{P}{The estimated joint excursion probability of each region.}
#' \item{P.err}{The Monte Carlo standard errors of `P`.}
#' \item{labels}{A vector with the region of each node, where 0 means that
#' the node is not in a region.}
#' \item{E}{The largest region, as an indicator vector.}
#' \item{F}{The excursion functions of the regions, as a sparse matrix with
#' one column for each region. See the details.}
#' \item{rho}{Marginal excursion probabilities.}
#' \item{mean}{The mean `mu`.}
#' \item{vars}{Marginal variances.}
#' \item{meta}{A list containing various information about the calculation.}
#' @export
#' @details
#' A connected region can only contain nodes with marginal excursion
#' probability at least \eqn{1-\alpha}{1-alpha}, so the search is done
#' within each connected component of these nodes. In each component, a
#' region is grown from a start point, which gives a sequence of connected
#' regions where each region contains the previous one. The joint excursion
#' probabilities of all regions in the sequence are computed with one run of
#' the sequential importance sampler, and the largest region with
#' probability at least \eqn{1-\alpha}{1-alpha} is kept. The start points
#' are the local maxima of the marginal excursion probabilities in the
#' component, computed on the standardised scale
#' \eqn{(\mu - u)/\sigma}{(mu - u)/sigma} to avoid ties for probabilities
#' close to one. Up to `n.starts` local maxima are tried, largest first,
#' since the largest maximum does not always give the largest region, and
#' the largest region over the start points is kept. Local maxima in the
#' largest region found so far are skipped. Finding the largest connected
#' region is a hard combinatorial problem, and the region that is found is
#' not guaranteed to be the largest one. The joint probability of each returned
#' region is however computed as in [excursions()].
#'
#' With `growth = 'bound'`, the node added in each step is the one that
#' decreases the Hunter lower bound of the joint probability the least. This
#' uses the pairwise failure probabilities of neighbouring nodes, which
#' are computed from the marginal variances and the covariances between
#' neighbours. These are obtained from the selected inverse of `Q`. Edges of
#' `graph` that are not in the sparsity pattern of the Cholesky factor of
#' `Q` are used for connectivity, but not for the bound.
#'
#' Each region has an excursion function, which is given for the nodes of
#' the region and for their neighbours. A region \eqn{R}{R} is the largest
#' set in a sequence of growing connected sets, and the excursion function
#' at a node of the region is the joint excursion probability of the
#' smallest set in the sequence that contains the node, so the region is the
#' set where the excursion function is at least \eqn{1-\alpha}{1-alpha}.
#' At a neighbour \eqn{j}{j} of the region, the excursion function is the
#' joint excursion probability of \eqn{R}{R} together with \eqn{j}{j}. These
#' are computed for all neighbours with one run of the sampler, and if a
#' neighbour that could be in the region has a probability of at least
#' \eqn{1-\alpha}{1-alpha}, it is added to the region. The excursion
#' function is zero at the nodes of the other regions and at all other
#' nodes, and it is used by [continuous()] to compute continuous domain
#' regions.
#'
#' After a region is found, its nodes are removed, and the remaining nodes of
#' its component are split into new components that are searched in the same
#' way. Regions grown from the other start points of the component
#' are reused in the new components if they do not overlap the removed
#' region, so these start points do not need to be computed again. Every node 
#' with marginal excursion probability at least \eqn{1-\alpha}{1-alpha} is a 
#' region of size one, so without `max.regions` or `min.size`, the regions 
#' cover all such nodes.
#'
#' The regions are computed jointly for all nodes in `mu`, with the other
#' nodes integrated out, so each start point requires a Cholesky
#' factorisation of `Q`. Components with a single node are computed directly from the
#' marginal probabilities.
#' @author David Bolin \email{davidbolin@@gmail.com}
#' @references Bolin, D. and Lindgren, F. (2015) *Excursion and contour
#' uncertainty regions for latent Gaussian models*, JRSS-series B, vol 77,
#' no 1, pp 85-106.
#' @seealso [excursions()]
#'
#' @examples
#' ## A field on a line with two bumps
#' n <- 50
#' Q <- sparseMatrix(
#'   i = c(1:n, 2:n), j = c(1:n, 1:(n - 1)),
#'   x = c(1, rep(1 + 0.9^2, n - 2), 1, rep(-0.9, n - 1)) / (1 - 0.9^2),
#'   dims = c(n, n), symmetric = TRUE
#' )
#' x <- seq(0, 1, length.out = n)
#' mu <- 4 * exp(-(x - 0.3)^2 / 0.01) + 3 * exp(-(x - 0.75)^2 / 0.005) - 1
#' res <- excursions.regions(
#'   alpha = 0.1, u = 0, mu = mu, Q = Q, type = ">",
#'   min.size = 3, seed = 1, max.threads = 1
#' )
#' res$regions
#' res$P
#' plot(x, mu, type = "l")
#' points(x, mu, col = res$labels + 1, pch = 19)
excursions.regions <- function(alpha,
                               u,
                               mu,
                               Q,
                               type,
                               graph,
                               n.iter = 20000,
                               vars,
                               rho,
                               method = "EB",
                               ind,
                               max.regions = Inf,
                               min.size = 1,
                               growth = c("bound", "rho"),
                               n.starts = 10,
                               min.prominence = 0,
                               max.threads = 0,
                               seed,
                               verbose = 0,
                               tol = NULL,
                               size.tol = 0.001) {
  growth <- match.arg(growth)
  if (method == "QC") {
    qc <- TRUE
  } else if (method == "EB") {
    qc <- FALSE
  } else {
    stop("only EB and QC methods are supported.")
  }
  if (missing(alpha)) {
    stop("Must specify error probability")
  }
  if (missing(u)) {
    stop("Must specify level")
  }
  if (missing(mu)) {
    stop("Must specify mean value")
  }
  if (missing(Q)) {
    stop("Must specify a precision matrix")
  }
  if (missing(type)) {
    stop("Must specify type of excursion set")
  }
  if (!(type %in% c(">", "<"))) {
    stop("Only the types '>' and '<' are supported.")
  }
  if (qc && missing(rho)) {
    stop("rho must be provided if QC is used.")
  }
  if (missing(seed)) {
    seed <- NULL
  }

  mu <- private.as.vector(mu)
  n <- length(mu)
  Q <- private.as.dgCMatrix(Q)
  if (!all(dim(Q) == n)) {
    stop("The dimensions of Q do not match the length of mu.")
  }
  if (missing(graph)) {
    graph <- NULL
  }
  G <- private.regions.graph(graph, Q, n)
  if (missing(vars)) {
    vars <- NULL
  } else {
    vars <- private.as.vector(vars)
  }
  if (missing(rho)) {
    rho <- NULL
  } else {
    rho <- private.as.vector(rho)
  }
  if (missing(ind)) {
    ind <- NULL
  }

  out <- private.regions(
    alpha = alpha, u = u, configs = list(list(mu = mu, Q = Q, vars = vars)),
    weights = 1, type = type, qc = qc, rho = rho, G = G, ind = ind,
    n.iter = n.iter, max.regions = max.regions, min.size = min.size,
    growth = growth, n.starts = n.starts, min.prominence = min.prominence,
    max.threads = max.threads, seed = seed, verbose = verbose, tol = tol,
    size.tol = size.tol
  )
  out$meta <- list(
    calculation = "regions",
    type = type,
    level = u,
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
    ind = NULL,
    call = match.call()
  )
  class(out) <- "excurobj"
  out
}

## The search for connected excursion regions. The configurations are
## lists with the mean mu, the precision Q and optionally the variances vars,
## and the joint probabilities are mixtures over the configurations with the
## given weights. The growth of the regions uses the marginal probabilities
## of the mixture and the correlations of the first configuration.
private.regions <- function(alpha, u, configs, weights, type, qc, rho, G, ind,
                            n.iter, max.regions, min.size, growth, n.starts,
                            min.prominence, max.threads, seed, verbose,
                            tol = NULL, size.tol = NULL) {
  n <- length(configs[[1]]$mu)
  n.conf <- length(configs)
  weights <- weights / sum(weights)

  ## Selected inverse of the first configuration. The edges of the graph are
  ## added to its pattern, so that their covariances are computed.
  if (growth == "bound" || is.null(configs[[1]]$vars)) {
    if (verbose) {
      cat("Calculate selected inverse\n")
    }
    sel <- private.selected.inverse(configs[[1]]$Q, G, max.threads = max.threads)
    if (is.null(configs[[1]]$vars)) {
      configs[[1]]$vars <- sel$vars
    }
  }

  ## Limits of each configuration. On the standardised scale, node i fails
  ## to exceed the level if Z_i < h_i.
  for (j in seq_len(n.conf)) {
    cf <- configs[[j]]
    if (is.null(cf$vars)) {
      cf$vars <- excursions.variances(Q = cf$Q)
    }
    if (is.null(rho)) {
      cf$marg <- excursions.marginals(
        type = type, vars = cf$vars, mu = cf$mu, u = u, QC = qc
      )
    } else {
      cf$marg <- excursions.marginals(
        type = type, rho = rho, vars = cf$vars, mu = cf$mu, u = u, QC = qc
      )
    }
    cf$limits <- excursions.setlimits(cf$marg, cf$vars, type, QC = qc, u, cf$mu)
    if (type == ">") {
      cf$h <- cf$limits$a / sqrt(cf$vars)
    } else {
      cf$h <- -cf$limits$b / sqrt(cf$vars)
    }
    configs[[j]] <- cf
  }

  ## The marginal failure probabilities, mixed over the configurations, and
  ## the corresponding levels h on the standardised scale
  if (n.conf == 1) {
    h <- configs[[1]]$h
    fail <- pnorm(h)
  } else {
    fail <- numeric(n)
    for (j in seq_len(n.conf)) {
      fail <- fail + weights[j] * pnorm(configs[[j]]$h)
    }
    h <- qnorm(fail)
  }
  cand <- fail <= alpha
  if (!is.null(ind)) {
    cand <- cand & private.ind.logical(ind, n)
  }
  sd <- sqrt(configs[[1]]$vars)

  ## Pairwise failure probabilities of the edges, in the order of G@i
  G.row <- G@i + 1L
  G.col <- rep(seq_len(n), diff(G@p))
  q <- numeric(length(G.row))
  if (growth == "bound") {
    e <- which(cand[G.row] & cand[G.col] & G.row < G.col &
      is.finite(h[G.row]) & is.finite(h[G.col]))
    r <- sel$cov(G.row[e], G.col[e]) / (sd[G.row[e]] * sd[G.col[e]])
    known <- !is.na(r)
    qe <- numeric(length(e))
    qe[known] <- .Call("regions_bvn_lower",
      as.double(h[G.row[e][known]]), as.double(h[G.col[e][known]]),
      as.double(r[known]),
      PACKAGE = "excursions"
    )
    key.e <- G.row[e] * (n + 1) + G.col[e]
    key <- pmin(G.row, G.col) * (n + 1) + pmax(G.row, G.col)
    q <- qe[match(key, key.e)]
    q[is.na(q)] <- 0
  }

  components <- function(mask) {
    .Call("regions_components", G@p, G@i, as.integer(mask),
      PACKAGE = "excursions"
    )
  }

  ## Start points in a component: the local maxima of the standardised
  ## margin -h with prominence at least min.prominence, largest first. A node
  ## is a local maximum if no neighbour in the component is larger, and ties
  ## are broken by the node index, so that a plateau gives a single start
  ## point. The largest maximum has infinite prominence, so it is always a
  ## start point.
  starts <- function(nodes, mask) {
    z <- -h
    e <- mask[G.row] & mask[G.col]
    r <- G.row[e]
    s <- G.col[e]
    beaten <- z[s] > z[r] | (z[s] == z[r] & s < r)
    lmax <- setdiff(nodes, r[beaten])
    if (min.prominence > 0 && length(lmax) > 1) {
      prom <- .Call("regions_prominence", G@p, G@i, as.double(z),
        as.integer(nodes[order(-z[nodes], nodes)]),
        PACKAGE = "excursions"
      )
      prom[is.nan(prom)] <- 0
      lmax <- lmax[prom[lmax] >= min.prominence]
    }
    lmax[order(-z[lmax])]
  }

  ## The largest region in a component that is grown from a start point. The
  ## nodes are ordered by the growth order, and the region is the largest
  ## prefix of that order with joint probability at least 1 - alpha. The
  ## result also contains the nodes that determine it, which are the region
  ## and the next node in the growth order.
  grow <- function(start, mask) {
    ord <- .Call("regions_grow", G@p, G@i, as.double(q), as.double(fail),
      as.double(-h), as.integer(mask), as.integer(start), growth == "bound",
      PACKAGE = "excursions"
    )
    K <- length(ord)
    ## The first node of the growth order is integrated first, that is, it is
    ## placed last. The configurations have the same sparsity pattern.
    cind <- integer(n)
    cind[ord] <- rev(seq_len(K))
    reo <- private.camd(configs[[1]]$Q, cind)
    Pk <- E2 <- numeric(K)
    for (j in seq_len(n.conf)) {
      a <- rep(-Inf, n)
      b <- rep(Inf, n)
      a[ord] <- configs[[j]]$limits$a[ord]
      b[ord] <- configs[[j]]$limits$b[ord]
      ## The mixture can only be at least 1 - alpha while the probability of
      ## this configuration is at least 1 - alpha / weight, so the
      ## integration can stop there
      res <- excursions.call(a, b, reo, configs[[j]]$Q,
        lim = max(0, 1 - alpha / weights[j]), K = n.iter, max.size = K,
        n.threads = max.threads, seed = seed,
        tol = tol, tol.level = 1 - alpha, size.tol = size.tol
      )
      Pk <- Pk + weights[j] * res$Pv[n - seq_len(K) + 1]
      E2 <- E2 + (weights[j] * res$Ev[n - seq_len(K) + 1])^2
    }
    Ek <- sqrt(E2)
    k <- sum(Pk >= 1 - alpha)
    if (k == 0) {
      ## Only possible due to Monte Carlo error, since the first node has
      ## marginal probability at least 1 - alpha
      k <- 1
    }
    list(
      start = start, region = ord[seq_len(k)], used = ord[seq_len(min(k + 1, K))],
      F = Pk[seq_len(k)], P = Pk[k], P.err = Ek[k]
    )
  }

  is.better <- function(res, best) {
    is.null(best) || length(res$region) > length(best$region) ||
      (length(res$region) == length(best$region) && res$P > best$P)
  }

  ## The largest region in a component over the start points, where ties
  ## are broken by the joint probability. The runs are earlier results for
  ## start points in the component that are still valid, and count as tried
  ## start points. Start points in the largest region found so far are
  ## skipped, since growing from them mostly retraces that region, and they
  ## do not count towards n.starts. Returns the largest region and all runs.
  evaluate <- function(nodes, runs) {
    if (length(nodes) == 1) {
      return(list(
        best = list(
          region = nodes, F = 1 - fail[nodes], P = 1 - fail[nodes], P.err = 0
        ),
        runs = list()
      ))
    }
    mask <- logical(n)
    mask[nodes] <- TRUE
    best <- NULL
    for (res in runs) {
      if (is.better(res, best)) {
        best <- res
      }
    }
    tried <- vapply(runs, function(x) x$start, 0L)
    for (start in starts(nodes, mask)) {
      if (length(runs) >= n.starts ||
        (!is.null(best) && length(best$region) == length(nodes))) {
        break
      }
      if (start %in% tried || (!is.null(best) && start %in% best$region)) {
        next
      }
      res <- grow(start, mask)
      runs[[length(runs) + 1]] <- res
      if (is.better(res, best)) {
        best <- res
      }
    }
    list(best = best, runs = runs)
  }

  ## The probabilities P(R and node j) for a region R, ordered by its growth
  ## order, and the nodes j in layer, with standard errors. The layer nodes
  ## are sampled without constraint after the region, see shapeIntProbe.
  probe <- function(region, layer) {
    K <- length(region)
    cind <- integer(n)
    cind[layer] <- 1L
    cind[region] <- 1L + rev(seq_len(K))
    reo <- private.camd(configs[[1]]$Q, cind)
    is.probe <- integer(n)
    is.probe[layer] <- 1L
    Pl <- E2 <- numeric(length(layer))
    for (j in seq_len(n.conf)) {
      a <- rep(-Inf, n)
      b <- rep(Inf, n)
      a[region] <- configs[[j]]$limits$a[region]
      b[region] <- configs[[j]]$limits$b[region]
      res <- private.regions.probe.call(a, b,
        pa = configs[[j]]$limits$a, pb = configs[[j]]$limits$b,
        probe = is.probe, reo = reo, Q = configs[[j]]$Q, K = n.iter,
        max.size = K + length(layer), n.threads = max.threads, seed = seed,
        tol = tol
      )
      Pl <- Pl + weights[j] * res$P[layer]
      E2 <- E2 + (weights[j] * res$E[layer])^2
    }
    list(P = Pl, P.err = sqrt(E2))
  }

  ## The neighbours of a set of nodes that can be in a region
  allowed <- if (is.null(ind)) rep(TRUE, n) else private.ind.logical(ind, n)
  neighbours <- function(nodes) {
    mask <- logical(n)
    mask[nodes] <- TRUE
    nb <- unique(G.row[mask[G.col]])
    nb[!mask[nb] & allowed[nb]]
  }

  ## Final region of a component, and its excursion function one layer
  ## outside the region. The excursion function at a node in the region is
  ## the joint probability of the region up to that node in the growth order,
  ## and at a neighbour j of the region it is P(R and node j). If a
  ## neighbour in the component has P(R and node j) >= 1 - alpha, the node
  ## with the largest probability is added to the region, and the layer is
  ## recomputed.
  finalize <- function(res, nodes) {
    region <- res$region
    Fr <- res$F
    P <- res$P
    P.err <- res$P.err
    repeat {
      layer <- neighbours(region)
      if (length(layer) == 0) {
        Fl <- numeric(0)
        break
      }
      pr <- probe(region, layer)
      Fl <- pr$P
      ext <- which(layer %in% nodes & Fl >= 1 - alpha)
      if (length(ext) == 0) {
        break
      }
      i <- ext[which.max(Fl[ext])]
      ## The probability cannot increase when a node is added, but the
      ## estimates come from different runs of the sampler
      P <- min(Fl[i], P)
      P.err <- pr$P.err[i]
      region <- c(region, layer[i])
      Fr <- c(Fr, P)
    }
    list(
      region = region, F = Fr, layer = layer, F.layer = pmin(Fl, P),
      P = P, P.err = P.err
    )
  }

  ## Best-first search over the components. An unevaluated component has the
  ## upper bound given by its size, and an evaluated component has the size
  ## of its region. Removing a region splits the rest of its component into
  ## new components, while the other components are unchanged.
  ##
  ## The runs from the other start points of the component are kept if they
  ## do not use any node of the removed region. Growing from the same start
  ## point in the new component then gives the same growth order up to the
  ## node after the region, so the run gives the same region, and its joint
  ## probability does not depend on the removed nodes.
  pending <- list()
  add.components <- function(mask, runs = list()) {
    lab <- components(mask)
    run.lab <- lab[vapply(runs, function(x) x$start, 0L)]
    for (l in unique(lab[lab > 0])) {
      nodes <- which(lab == l)
      if (length(nodes) >= min.size) {
        pending[[length(pending) + 1]] <<- list(
          nodes = nodes, size = length(nodes), result = NULL,
          runs = runs[run.lab == l]
        )
      }
    }
  }
  add.components(cand)

  regions <- Fs <- list()
  P <- P.err <- numeric(0)
  while (length(regions) < max.regions && length(pending) > 0) {
    size <- vapply(pending, function(x) x$size, 0)
    done <- vapply(pending, function(x) !is.null(x$result), TRUE)
    ## Among the components with the largest bound, prefer evaluated ones
    i <- order(-size, !done)[1]
    if (!done[i]) {
      if (verbose) {
        cat("Evaluate component with", size[i], "nodes\n")
      }
      ev <- evaluate(pending[[i]]$nodes, pending[[i]]$runs)
      pending[[i]]$result <- ev$best
      pending[[i]]$runs <- ev$runs
      pending[[i]]$size <- length(ev$best$region)
      if (pending[[i]]$size < min.size) {
        pending[[i]] <- NULL
      }
      next
    }
    comp <- pending[[i]]
    pending[[i]] <- NULL
    fin <- finalize(comp$result, comp$nodes)
    regions[[length(regions) + 1]] <- fin$region
    Fs[[length(Fs) + 1]] <- list(
      nodes = c(fin$region, fin$layer), F = c(fin$F, fin$F.layer)
    )
    P <- c(P, fin$P)
    P.err <- c(P.err, fin$P.err)
    rest <- setdiff(comp$nodes, fin$region)
    if (length(rest) >= min.size) {
      mask <- logical(n)
      mask[rest] <- TRUE
      valid <- vapply(comp$runs, function(x) all(mask[x$used]), TRUE)
      add.components(mask, comp$runs[valid])
    }
  }

  ## The regions are found in order of decreasing size, except that a
  ## region in the rest of a component can be larger than the region that
  ## was removed from it
  o <- order(-lengths(regions))
  regions <- lapply(regions[o], sort)
  P <- P[o]
  P.err <- P.err[o]
  Fs <- Fs[o]
  F_ <- sparseMatrix(
    i = as.integer(unlist(lapply(Fs, function(x) x$nodes))),
    j = rep(seq_along(Fs), vapply(Fs, function(x) length(x$nodes), 0L)),
    x = as.double(unlist(lapply(Fs, function(x) x$F))),
    dims = c(n, length(Fs))
  )
  labels <- integer(n)
  for (k in seq_along(regions)) {
    labels[regions[[k]]] <- k
  }
  E <- as.numeric(labels == 1)

  ## The excursion function of a region is zero at the nodes of the other
  ## regions, so that the region is the set where it is at least 1 - alpha
  if (length(regions) > 0) {
    tr <- summary(F_)
    other <- labels[tr$i] > 0 & labels[tr$i] != tr$j
    F_ <- sparseMatrix(
      i = tr$i[!other], j = tr$j[!other], x = tr$x[!other],
      dims = dim(F_)
    )
  }

  list(
    regions = regions,
    P = P,
    P.err = P.err,
    labels = labels,
    E = E,
    F = F_,
    rho = configs[[1]]$marg$rho,
    mean = configs[[1]]$mu,
    vars = configs[[1]]$vars
  )
}

## Symmetric neighbourhood graph as a dgCMatrix pattern without diagonal
private.regions.graph <- function(graph, Q, n) {
  if (is.null(graph)) {
    graph <- Q
  } else if (inherits(graph, c("fm_mesh_2d", "inla.mesh"))) {
    graph <- graph$graph$vv
  }
  if (!all(dim(graph) == n)) {
    stop("The dimensions of graph do not match the length of mu.")
  }
  tr <- private.sparse.gettriplet(graph)
  keep <- tr$x != 0 & tr$i != tr$j
  i <- tr$i[keep]
  j <- tr$j[keep]
  private.as.dgCMatrix(sparseMatrix(
    i = c(i, j), j = c(j, i), x = 1,
    dims = c(n, n)
  ))
}

## Selected inverse of Q: the marginal variances, and a function that returns
## the covariances of pairs of nodes, which is NA for pairs that are not in
## the sparsity pattern of the Cholesky factor. If the graph G is given, its
## edges are added to the pattern of Q as explicit zeros, which puts them in
## the pattern of the Cholesky factor without changing Q.
private.selected.inverse <- function(Q, G = NULL, max.threads = 0) {
  Q <- private.as.dgCMatrix(Q)
  if (!is.null(G)) {
    Qt <- as(Q, "TsparseMatrix")
    Gt <- as(private.as.dgCMatrix(G), "TsparseMatrix")
    Q <- as(methods::new("dgTMatrix",
      i = c(Qt@i, Gt@i), j = c(Qt@j, Gt@j),
      x = c(Qt@x, numeric(length(Gt@i))), Dim = dim(Q)
    ), "CsparseMatrix")
  }
  ## CHOLMOD chooses the supernodal factorization if it is faster. Its
  ## explicit zeros are kept, since they are part of the closed pattern that
  ## the recursion for the selected inverse needs.
  ch <- Matrix::Cholesky(Q, LDL = FALSE, perm = TRUE, super = NA)
  private.selected.inverse.factor(private.factor.lower(ch), ch@perm + 1L,
    max.threads = max.threads
  )
}

## Selected inverse from a Cholesky factor of Q, see
## private.selected.inverse. L is the lower triangular factor of Q[perm, perm],
## or of Q if perm is NULL. An upper triangular factor is transposed. The
## recursion uses max.threads threads, see excursions.variances.
private.selected.inverse.factor <- function(L, perm = NULL, max.threads = 0) {
  L <- private.as.dtCMatrix(L)
  if (L@uplo == "U") {
    L <- t(L)
  }
  if (L@diag == "U") {
    L <- as(L, "generalMatrix")
  }
  n <- nrow(L)
  if (is.null(perm)) {
    perm <- seq_len(n)
  }
  z <- .Call("Qinv_selected", L@p, L@i, as.double(L@x),
    as.integer(max.threads),
    PACKAGE = "excursions"
  )
  vars <- numeric(n)
  vars[perm] <- z[L@p[-(n + 1)] + 1L]
  iperm <- integer(n)
  iperm[perm] <- seq_len(n)
  ## L@i are the rows of the lower triangle, so each key is (column, row)
  ## with column <= row
  key <- rep(seq_len(n), diff(L@p)) * (n + 1) + (L@i + 1L)
  cov <- function(i, j) {
    a <- iperm[i]
    b <- iperm[j]
    z[match(pmin(a, b) * (n + 1) + pmax(a, b), key)]
  }
  list(vars = vars, cov = cov)
}

## Sequential importance sampling as in excursions.call, where the nodes with
## probe != 0 are probes, see shapeIntProbe. Returns the probabilities and
## standard errors of the probes in the original order.
private.regions.probe.call <- function(a, b, pa, pb, probe, reo, Q, K,
                                       max.size, n.threads, seed,
                                       tol = NULL) {
  n <- length(a)
  L <- suppressWarnings(private.Cholesky(Q[reo, reo], perm = FALSE)$R)
  finite <- function(x) {
    x <- x[reo]
    x[x == Inf] <- .Machine$double.xmax
    x[x == -Inf] <- -.Machine$double.xmax
    x
  }
  a <- finite(a)
  b <- finite(b)
  pa <- finite(pa)
  pb <- finite(pb)
  if (!is.null(seed)) {
    seed <- as.integer(seed)
    if (length(seed) == 1) {
      seed <- rep(seed, 6)
    }
    seed.provided <- 1L
  } else {
    seed <- integer(6)
    seed.provided <- 0L
  }
  L_ipx <- private.sparse.get_ipx(L)
  out <- .C("shapeIntProbe",
    Mp = as.integer(L_ipx$p), Mi = as.integer(L_ipx$i),
    Mv = as.double(L_ipx$x), a = as.double(a), b = as.double(b),
    opts = as.integer(c(n, K, max.size, n.threads, seed.provided)),
    lim = as.double(0), Pv = double(n), Ev = double(n), seed_in = seed,
    probe = as.integer(probe[reo]), pa = as.double(pa),
    pb = as.double(pb), Pp = double(n), Pe = double(n),
    ## As in gaussint, tol = 0 gives a single batch of K samples, and the
    ## first batch has 1000 samples otherwise. The error of the probability
    ## of the region is controlled, which bounds the errors of the probes.
    adapt = as.double(c(if (is.null(tol)) 0 else tol, 0, 1000, 0)),
    K_used = integer(1),
    PACKAGE = "excursions"
  )
  P <- E <- numeric(n)
  P[reo] <- out$Pp
  E[reo] <- out$Pe
  list(P = P, E = E)
}
