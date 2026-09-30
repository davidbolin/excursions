## utils.R
##
##   Copyright (C) 2013,2016,2020 David Bolin, Finn Lindgren
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


## Calculate upper triangular Cholesky decomposition, optionally with
## permutation. All Matrix::Cholesky options are allowed.
## Returns list(R=dtCMatrix, reo=integer vector, ireo=integer vector)
private.Cholesky <- function(A, ...) {
  L <- expand(Matrix::Cholesky(private.as.dgCMatrix(A), ...))
  n <- nrow(A)
  ireo <- integer(n)
  ireo[L$P@perm] <- seq_len(n)
  reo <- integer(n)
  reo[ireo] <- seq_len(n)
  list(R = private.as.dtCMatrixU(L$L), reo = reo, ireo = ireo)
}



#' Calculate variances from a sparse precision matrix
#'
#' `excursions.variances` calculates the diagonal of the inverse of a sparse
#' symmetric positive definite matrix `Q`.
#'
#' @param L Cholesky factor of precision matrix.
#' @param Q Precision matrix.
#' @param max.threads Not used. The computation is sequential, and the argument is
#' kept for backwards compatibility.
#'
#' @return A vector with the variances.
#' @export
#' @details The method for calculating the
#' diagonal requires the Cholesky factor, `L`, of `Q`, which should be supplied if
#' available. If `Q` is provided, the cholesky factor is
#' calculated and the variances are then returned in the same ordering as `Q`.
#' If `L` is provided, the variances are returned in the same ordering as `L`,
#' even if `L@invpivot` exists.
#' @author David Bolin \email{davidbolin@@gmail.com}
#'
#' @examples
#' ## Create a tridiagonal precision matrix
#' n <- 21
#' Q <- Matrix(toeplitz(c(1, -0.1, rep(0, n - 2))))
#' v2 <- excursions.variances(Q = Q, max.threads = 2)
#' ## var2 should be the same as:
#' v1 <- diag(solve(Q))
excursions.variances <- function(L, Q, max.threads = 0) {
  if (!missing(L) && !is.null(L)) {
    perm <- NULL
    if (inherits(L, "sparseMatrix") && !inherits(L, "triangularMatrix")) {
      ## For example a supernodal factor, which may store zeros on the other
      ## side of the diagonal
      L <- Matrix::drop0(L)
    }
    L <- private.as.dtCMatrix(L)
  } else {
    ## Keep the factor lower triangular, as the C code works on columns of L.
    ch <- Matrix::Cholesky(private.as.dgCMatrix(Q), LDL = FALSE, perm = TRUE)
    perm <- ch@perm + 1L
    L <- as(ch, "CsparseMatrix")
  }
  lower <- L@uplo == "L"
  if (L@diag == "U") {
    L <- as(L, "generalMatrix")
  }

  variances <- .Call("Qinv", L@p, L@i, as.double(L@x), lower)

  if (is.null(perm)) {
    variances
  } else {
    variances[perm] <- variances
    variances
  }
}


excursions.marginals <- function(type, rho, vars, mu, u, QC = FALSE) {
  rl <- list()
  if (type == "=" || type == "!=") {
    if (QC) {
      rl$rho_ngu <- rho
      rl$rho_l <- pnorm(mu - u, sd = sqrt(vars), lower.tail = FALSE)
      rl$rho_u <- 1 - rl$rho_l
      rl$rho <- pmax(rl$rho_u, rl$rho_l)
      rl$rho_ng <- pmax(rl$rho_ngu, 1 - rl$rho_ngu)
    } else {
      if (!missing(rho)) {
        rl$rho_u <- rho
      } else {
        rl$rho_u <- 1 - pnorm(mu - u, sd = sqrt(vars), lower.tail = FALSE)
      }
      rl$rho_l <- 1 - rl$rho_u
      rl$rho <- pmax(rl$rho_u, rl$rho_l)
    }
  } else {
    if (QC) {
      rl$rho_ng <- rho
      if (type == ">") {
        rl$rho <- pnorm(mu - u, sd = sqrt(vars))
      } else {
        rl$rho <- pnorm(mu - u, sd = sqrt(vars), lower.tail = FALSE)
      }
    } else {
      if (missing(rho)) {
        if (type == ">") {
          rl$rho <- pnorm(mu - u, sd = sqrt(vars))
        } else {
          rl$rho <- pnorm(mu - u, sd = sqrt(vars), lower.tail = FALSE)
        }
      } else {
        rl$rho <- rho
      }
    }
  }
  return(rl)
}


excursions.permutation <- function(rho, ind, use.camd = TRUE, alpha, Q) {
  if (!missing(ind) && !is.null(ind)) {
    rho[!private.ind.logical(ind, length(rho))] <- -1
  }
  n <- length(rho)
  v.s <- sort(rho, index.return = TRUE)
  reo <- v.s$ix
  rho_sort <- v.s$x
  ireo <- integer(n)
  ireo[reo] <- 1:n
  if (use.camd) {
    k <- 0
    i <- n
    # add nodes to lower bound
    cindr <- cind <- rep(0, n)
    while (i > 0 && rho_sort[i] > 1 - alpha) {
      cindr[i] <- k
      i <- i - 1
      k <- k + 1
    }
    if (i > 0) {
      # reorder nodes below the lower bound for sparsity
      while (i > 0) {
        cindr[i] <- k
        i <- i - 1
      }
      # change back to original ordering
      cind <- k - cindr[ireo]
      reo <- private.camd(Q, cind)
    }
  }
  return(reo)
}

## Constrained approximate minimum degree ordering of Q, where the nodes with
## constraint cind == c are ordered before the nodes with cind == c + 1.
##
## A constraint set with a single node fixes the position of that node, and
## CAMD is slow when there are many constraint sets. Runs of consecutive
## single node sets are therefore merged into one set before calling CAMD,
## and the nodes of each merged run are then put back in their fixed order.
## This does not change the ordering, since the order of the nodes within a
## run does not change the graph that remains for the other sets.
private.camd <- function(Q, cind) {
  n <- length(cind)
  Q <- private.as.dgCMatrix(Q)
  Q_ipx <- private.sparse.get_ipx(Q)

  sets <- sort(unique(cind))
  single <- tabulate(match(cind, sets), length(sets)) == 1
  ## New set index, shared by consecutive single node sets
  merged <- cumsum(!(single & c(FALSE, single[-length(single)])))
  cind.merged <- merged[match(cind, sets)] - 1L

  out <- .C("reordering",
    nin = as.integer(n), Mp = as.integer(Q_ipx$p),
    Mi = as.integer(Q_ipx$i), reo = integer(n),
    cind = as.integer(cind.merged)
  )
  reo <- out$reo + 1L

  ## Restore the fixed order within each merged run. The nodes of a set are
  ## consecutive in reo, so sorting by the original set index within the
  ## positions of the merged set restores it.
  runs <- which(tabulate(cind.merged + 1L, max(merged)) > 1 &
    vapply(split(single, merged), all, TRUE))
  for (r in runs) {
    pos <- which(cind.merged[reo] == r - 1L)
    reo[pos] <- reo[pos][order(cind[reo[pos]])]
  }
  reo
}


excursions.setlimits <- function(marg, vars, type, QC, u, mu) {
  if (QC) {
    if (type == "<") {
      uv <- sqrt(vars) * qnorm(pmin(pmax(marg$rho_ng, 0), 1))
    } else if (type == ">") {
      uv <- sqrt(vars) * qnorm(pmin(pmax(marg$rho_ng, 0), 1), lower.tail = FALSE)
    } else if (type == "=" || type == "!=") {
      uv <- sqrt(vars) * qnorm(pmin(pmax(marg$rho_ngu, 0), 1), lower.tail = FALSE)
    }
  } else {
    uv <- u - mu
  }
  if (type == "=" || type == "!=") {
    if (QC) {
      a <- b <- uv
      a[marg$rho_ngu <= 0.5] <- -Inf
      b[marg$rho_ngu > 0.5] <- Inf
    } else {
      a <- b <- uv
      a[marg$rho_u <= 0.5] <- -Inf
      b[marg$rho_u > 0.5] <- Inf
    }
  } else if (type == ">") {
    a <- uv
    b <- rep(Inf, length(mu))
  } else if (type == "<") {
    a <- rep(-Inf, length(mu))
    b <- uv
  }

  list(a = a, b = b)
}



excursions.call <- function(a, b, reo, Q, is.chol = FALSE, lim, K, max.size, n.threads, seed) {
  if (is.chol && !identical(as.integer(reo), seq_len(length(reo)))) {
    ## The factor is for the original ordering, so form the precision matrix
    ## and factorise it in the integration order
    L <- private.as.dtCMatrixU(Q)
    Q <- crossprod(L)
    is.chol <- FALSE
  }
  if (!is.chol) {
    a.sort <- a[reo]
    b.sort <- b[reo]
    Q <- Q[reo, reo]

    L <- suppressWarnings(private.Cholesky(Q, perm = FALSE)$R)

    res <- gaussint(
      Q.chol = L, a = a.sort, b = b.sort, lim = lim,
      n.iter = K, max.size = max.size,
      max.threads = n.threads, seed = seed
    )
  } else {
    ## The integration order is the order of the factor
    res <- gaussint(
      Q.chol = Q, a = a, b = b, lim = lim, n.iter = K,
      max.size = max.size,
      max.threads = n.threads, seed = seed
    )
  }
  return(res)
}


private.check.integer <- function(v) {
  if (is.null(v)) {
    stop("Anticipated scalar value, got NULL")
  } else if (!is.null(dim(v))) {
    stop("Anticipated scalar value, got matrix")
  } else if (length(v) > 1) {
    stop("Anticipated scalar value, got vector")
  }
}

## Logical vector of length n from indices, which may be given as a logical
## vector or as integer indices
private.ind.logical <- function(ind, n) {
  if (is.logical(ind)) {
    return(rep_len(ind, n))
  }
  lind <- rep(FALSE, n)
  lind[ind] <- TRUE
  lind
}

private.as.vector <- function(v) {
  if (is.null(v) || is.vector(v)) {
    return(v)
  }
  if (min(dim(v) > 1)) {
    stop("vector has wrong dimensions")
  }
  as.vector(v)
}

private.sparse.gettriplet <- function(M) {
  ## Get unique triplet representation:
  M <- private.as.dgTMatrix(M)
  ## Extract triplets:
  list(i = M@i + 1L, j = M@j + 1L, x = M@x)
}

private.sparse.get_ipx <- function(M) {
  if (!inherits(M, "CsparseMatrix")) {
    stop("M must be a CsparseMatrix object")
  }
  M <- private.as.dgCMatrix(M)
  ## Extract i,p,x in 0-based format:
  ## If M is a unit diagonal matrix, may have length(i)==0
  if (inherits(M, "triangularMatrix") &&
    (M@diag == "U") &&
    (length(M@i) == 0)) {
    list(
      i = seq_len(nrow(M)) - 1L,
      p = seq_len(nrow(M) + 1) - 1L,
      x = rep(1.0, nrow(M))
    )
  } else {
    list(i = M@i, p = M@p, x = M@x)
  }
}

private.as.dgTMatrix <- function(M, make_unique = TRUE) {
  if (is.null(M) || (!make_unique && is(M, "dgTMatrix"))) {
    return(M)
  } else {
    ## Convert into dgTMatrix format of Matrix. Make sure the
    ## representation is unique (ie no double triplets etc)
    ## convert through the 'dgCMatrix'-class to make it unique;
    return(as(private.as.dgCMatrix(M), "TsparseMatrix"))
  }
}

private.as.dgCMatrix <- function(M) {
  if (is.null(M) || is(M, "dgCMatrix")) {
    return(M)
  }

  if (!inherits(M, "Matrix")) {
    M <- as(M, "Matrix")
  }
  ## Convert into dgCMatrix format of Matrix.
  ## Convert via virtual class CsparseMatrix;
  ## this allows more general conversions than direct conversion.
  as(as(as(M, "dMatrix"), "generalMatrix"), "CsparseMatrix")
}

private.as.dtCMatrix <- function(M) {
  if (is.null(M) || is(M, "dtCMatrix")) {
    return(M)
  }

  if (!inherits(M, "Matrix")) {
    M <- as(M, "Matrix")
  }
  ## Convert into dtCMatrix format of Matrix.
  ## Convert via virtual class CsparseMatrix;
  ## this allows more general conversions than direct conversion.
  as(as(as(M, "dMatrix"), "triangularMatrix"), "CsparseMatrix")
}


## Transpose a lower triangular matrix into upper triangular dtCMatrix
## If already upper triangular dtCMatrix, the matrix is returned unchanged
private.as.dtCMatrixU <- function(M) {
  M <- private.as.dtCMatrix(M)
  if (M@uplo == "L") {
    t(M)
  } else {
    M
  }
}



##
# Quantile function of Gaussian mixture
##
##
# Quantile function of Gaussian mixtures, vectorised over locations.
# mu and sd are K x n matrices, with one column per location, and w are the
# K mixture weights. Returns the p-quantile at each of the n locations.
#
# The quantile lies between the smallest and largest component quantiles,
# which gives a starting bracket. Safeguarded Newton steps are then taken
# for all locations at once. Quantiles outside br are set to the nearest end point of br.
##
Fmix_inv_vec <- function(p, mu, sd, w, br = c(-1000, 1000),
                         tol = 1e-10, max.iter = 100) {
  mu <- as.matrix(mu)
  sd <- as.matrix(sd)
  n <- ncol(mu)
  zp <- stats::qnorm(p)
  qk <- mu + sd * zp
  lo <- apply(qk, 2, min)
  hi <- apply(qk, 2, max)
  x <- (lo + hi) / 2
  ## Quantiles beyond an end point of br are set to that end point
  Fmix.at <- function(y, idx) {
    colSums(w * stats::pnorm((y - mu[, idx, drop = FALSE]) /
      sd[, idx, drop = FALSE]))
  }
  below.br <- above.br <- logical(n)
  idx <- which(lo < br[1])
  below.br[idx] <- Fmix.at(br[1], idx) >= p
  idx <- which(hi > br[2])
  above.br[idx] <- Fmix.at(br[2], idx) <= p
  lo <- pmax(lo, br[1])
  hi <- pmin(hi, br[2])
  x[below.br] <- br[1]
  x[above.br] <- br[2]
  active <- which(hi > lo & !below.br & !above.br)
  dx.old <- hi - lo
  iter <- 0
  while (length(active) > 0 && iter < max.iter) {
    iter <- iter + 1
    xa <- x[active]
    z <- (rep(xa, each = nrow(mu)) - mu[, active, drop = FALSE]) /
      sd[, active, drop = FALSE]
    Fx <- colSums(w * stats::pnorm(z)) - p
    fx <- colSums(w * stats::dnorm(z) / sd[, active, drop = FALSE])
    below <- Fx < 0
    lo[active[below]] <- xa[below]
    hi[active[!below]] <- xa[!below]
    ## Newton step, unless it leaves the bracket or does not reduce the
    ## step length fast enough, in which case bisect (as in rtsafe)
    dx <- Fx / fx
    dx[Fx == 0] <- 0
    converged <- Fx == 0 | abs(dx) <= tol * (1 + abs(xa))
    xn <- xa - dx
    bisect <- !converged &
      (!is.finite(xn) | xn <= lo[active] | xn >= hi[active] |
        abs(2 * Fx) > abs(dx.old[active] * fx))
    xn[bisect] <- (lo[active[bisect]] + hi[active[bisect]]) / 2
    dx[bisect] <- xa[bisect] - xn[bisect]
    dx.old[active] <- dx
    x[active] <- xn
    done <- converged | abs(dx) <= tol * (1 + abs(xa))
    active <- active[!done]
  }
  x
}

##
# Function for optimization of interval for mixtures
##
fmix.opt <- function(x,
                     alpha,
                     sd,
                     Q.chol,
                     w,
                     mu,
                     limits,
                     verbose,
                     max.threads,
                     ind,
                     n.iter = 10000,
                     seed = NULL) {
  K <- dim(mu)[1]
  q.a <- Fmix_inv_vec(x / 2, mu = mu, sd = sd, w = w, br = limits)
  q.b <- Fmix_inv_vec(1 - x / 2, mu = mu, sd = sd, w = w, br = limits)

  prob <- 0
  stopped <- 0

  k.seq <- sort(w, decreasing = TRUE, index.return = TRUE)$ix

  ki <- 1
  for (k in k.seq) {
    ws <- 0
    if (ki < K) {
      ws <- sum(w[k.seq[(ki + 1):K]])
    }
    ki <- ki + 1
    lim <- (1 - alpha - prob - ws) / w[k]
    p <- gaussint(
      mu = mu[k, ],
      Q.chol = Q.chol[[k]],
      a = q.a,
      b = q.b,
      ind = ind,
      lim = max(0, lim),
      n.iter = n.iter,
      max.threads = max.threads,
      seed = seed
    )
    if (p$P == 0) {
      stopped <- 1
      break
    } else {
      prob <- prob + w[k] * p$P
    }
  }


  if (stopped == 1) { # too large alpha
    if (prob == 0) {
      val <- 10 * (1 + x)
    } else {
      val <- 1 + (prob - (1 - alpha))^2
    }
  } else { # too small x
    val <- (prob - (1 - alpha))^2
  }

  if (verbose) {
    cat("in optimization: ", x, " ", prob, " ", val, "\n")
  }

  val
}


fmix.samp.opt <- function(x,
                          alpha,
                          mu,
                          sd,
                          w,
                          limits,
                          samples,
                          verbose = FALSE) {
  q.a <- Fmix_inv_vec(x / 2, mu = mu, sd = sd, w = w, br = limits)
  q.b <- Fmix_inv_vec(1 - x / 2, mu = mu, sd = sd, w = w, br = limits)

  ## samples has one row per sample, compare column-wise on the transpose
  cover <- colSums(t(samples) > q.b | t(samples) < q.a) == 0

  prob <- mean(cover)
  val <- (prob - (1 - alpha))^2
  if (verbose) {
    cat("in optimization: ", x, " ", prob, " ", val, "\n")
  }

  return(val)
}

## Type 7 quantiles, as in stats::quantile, of each row of a matrix whose
## rows are sorted and free of NA.
private.row.quantile <- function(sorted, p) {
  index <- 1 + max(ncol(sorted) - 1, 0) * p
  lo <- floor(index)
  hi <- ceiling(index)
  qs <- sorted[, lo]
  if (index > lo) {
    x.hi <- sorted[, hi]
    i <- which(x.hi != qs)
    h <- index - lo
    qs[i] <- (1 - h) * qs[i] + h * x.hi[i]
  }
  qs
}

## Row-wise quantiles of samples, using the row-sorted samples if available.
private.samples.quantile <- function(samples, p, sorted = NULL) {
  if (is.null(sorted)) {
    apply(samples, 1, quantile, 1, probs = p)
  } else {
    private.row.quantile(sorted, p)
  }
}

## Sort each row of samples, or NULL if samples has missing values.
private.sort.rows <- function(samples) {
  if (anyNA(samples)) {
    return(NULL)
  }
  sorted <- t(apply(samples, 1, sort))
  if (ncol(samples) == 1) {
    sorted <- t(sorted)
  }
  sorted
}

fsamp.opt <- function(x, samples, verbose = FALSE, sorted = NULL) {
  q.a <- private.samples.quantile(samples, x / 2, sorted)
  q.b <- private.samples.quantile(samples, 1 - x / 2, sorted)
  if (is.null(sorted)) {
    prob <- mean(apply((samples < q.b) * (samples > q.a), 2, prod))
  } else {
    ## As a double vector, so that mean() gives the same result as above
    prob <- mean(as.double(colSums(!((samples < q.b) & (samples > q.a))) == 0))
  }
  if (verbose) {
    cat("in optimization: ", x, " ", prob, "\n")
  }
  return(prob)
}


mix.sample <- function(n.samp = 1, mu, Q.chol, w) {
  K <- length(mu)
  n <- length(mu[[1]])
  idx <- sample(seq_len(K), n.samp, prob = w, replace = TRUE)
  idx <- sort(idx)
  n.idx <- numeric(K)
  n.idx[] <- 0
  for (i in 1:K) {
    n.idx[i] <- sum(idx == i)
  }
  samples <- c()
  for (k in 1:K) {
    if (n.idx[k] > 0) {
      xx <- mu[[k]] + solve(Q.chol[[k]], matrix(rnorm(n.idx[k] * n), n, n.idx[k]))
      samples <- rbind(samples, t(as.matrix(xx)))
    }
  }
  return(samples)
}

excursions.rand <- function(n, seed, n.threads = 1) {
  if (!missing(seed) && !is.null(seed)) {
    seed_provided <- 1
    seed.in <- seed
  } else {
    seed_provided <- 0
    seed.in <- as.integer(rep(0, 6))
  }

  x <- rep(0, n)
  opt <- c(n, n.threads, seed_provided)
  out <- .C("testRand",
    opt = as.integer(opt),
    x = as.double(x),
    seed_in = as.integer(seed.in)
  )

  return(out$x)
}



#' Warnings free loading of add-on packages
#'
#' Turn off all warnings for require(), to allow clean completion
#' of examples that require unavailable Suggested packages.
#'
#' @param package The name of a package, given as a character string.
#' @param lib.loc a character vector describing the location of R library trees
#' to search through, or `NULL`.  The default value of `NULL`
#' corresponds to all libraries currently known to `.libPaths()`.
#' Non-existent library trees are silently ignored.
#' @param character.only a logical indicating whether `package` can be
#' assumed to be a character string.
#'
#' @return `require.nowarnings` returns (invisibly) `TRUE` if it succeeds, otherwise `FALSE`
#' @details `require(package)` acts the same as
#' `require(package, quietly = TRUE)` but with warnings turned off.
#' In particular, no warning or error is given if the package is unavailable.
#' Most cases should use `requireNamespace(package, quietly = TRUE)` instead,
#' which doesn't produce warnings.
#' @seealso [require()]
#' @export
#' @examples
#' ## This should produce no output:
#' if (require.nowarnings(nonexistent)) {
#'   message("Package loaded successfully")
#' }
require.nowarnings <- function(package, lib.loc = NULL, character.only = FALSE) {
  if (!character.only) {
    package <- as.character(substitute(package))
  }
  suppressWarnings(
    require(package,
      lib.loc = lib.loc,
      quietly = TRUE,
      character.only = TRUE
    )
  )
}



excursions.marginals.mc <- function(X, type, rho, mu, u) {
  rl <- list()
  if (type == "=" || type == "!=") {
    if (!missing(rho)) {
      rl$rho_u <- rho
    } else {
      rl$rho_u <- 1 - rowMeans(X < u)
    }
    rl$rho_l <- 1 - rl$rho_u
    rl$rho <- pmax(rl$rho_u, rl$rho_l)
  } else {
    if (missing(rho)) {
      if (type == ">") {
        rl$rho <- rowMeans(X > u)
      } else {
        rl$rho <- rowMeans(X < u)
      }
    }
  }
  rl
}

mcint <- function(X,
                  a,
                  b,
                  ind) {
  if (missing(a)) {
    stop("Must specify lower integration limit")
  }

  if (missing(b)) {
    stop("Must specify upper integration limit")
  }

  n <- length(a)
  if (length(b) != n) {
    stop("Vectors with integration limits are of different length.")
  }

  if (!missing(ind) && !is.null(ind)) {
    ind <- private.ind.logical(ind, n)
    a[!ind] <- -Inf
    b[!ind] <- Inf
  }

  inside <- a < X & X < b
  if (anyNA(inside)) {
    Pv <- rowMeans(apply(apply(apply(inside, 2, rev), 2, cumprod), 2, rev))
  } else {
    ## Pv[i] is the fraction of samples inside the limits for all j >= i
    Pv <- numeric(n)
    alive <- rep(TRUE, ncol(inside))
    for (i in rev(seq_len(n))) {
      alive <- alive & inside[i, ]
      Pv[i] <- mean(alive)
    }
  }

  # Estimate of MC error, not implemented yet
  Ev <- rep(0, n)

  list(Pv = Pv, Ev = Ev, P = Pv[1], E = Ev[1])
}


#' Summarise excurobj objects
#'
#' Summary method for class "excurobj"
#'
#' @param object an object of class "excurobj", usually, a result of a call
#'   to [excursions()].
#' @param ... further arguments passed to or from other methods.
#' @export
#' @method summary excurobj
summary.excurobj <- function(object, ...) {
  out <- list()
  class(out) <- "summary.excurobj"
  out$calculation <- object$meta$calculation
  out$call <- object$meta$call
  if (object$meta$calculation == "excursions") {
    if (object$meta$type == ">") {
      out$computation <- "Positive excursion set, E_{u,alpha}^+"
    } else if (object$meta$type == "<") {
      out$computation <- "Negative excursion set, E_{u,alpha}^-"
    } else if (object$meta$type == "=") {
      out$computation <- "Contour credible region, E_{u,alpha}^c"
    } else {
      out$computation <- "Contour avoiding set, E_{u,\alpha}"
    }
    out$u <- object$meta$level
    out$alpha <- object$meta$alpha
    out$F.limit <- object$meta$F.limit
    out$method <- object$meta$method
  } else if (object$meta$calculation == "simconf") {
    out$computation <- "Simultaneous confidence band"
    out$alpha <- object$alpha
  } else if (object$meta$calculation == "contourmap") {
    out$computation <- "Contour map"
    out$u <- object$u
    out$type <- object$meta$contourmap.type
    out$F.computed <- object$meta$F.computed
    if (out$F.computed) {
      out$F.limit <- object$meta$F.limit
    }

    if (is.null(object$P0) && is.null(object$P1) && is.null(object$P2) &&
      is.null(object$P0.bound) && is.null(object$P1.bound) &&
      is.null(object$P2.bound)) {
    } else {
      out$measures <- list()
      if (!is.null(object$P0)) {
        out$measures$P0 <- object$P0
      }

      if (!is.null(object$P1)) {
        if (!is.null(object$P1.error)) {
          out$measures$P1 <- sprintf("%.4g (error %.5g)", object$P1, object$P1.error)
        } else {
          out$measures$P1 <- object$P1
        }
      }

      if (!is.null(object$P2)) {
        if (!is.null(object$P2.error)) {
          out$measures$P2 <- sprintf("%.4g (error %.5g)", object$P2, object$P2.error)
        } else {
          out$measures$P2 <- object$P2
        }
      }

      if (!is.null(object$P0.bound)) {
        out$measures$P0.bound <- object$P0.bound
      }

      if (!is.null(object$P1.bound)) {
        out$measures$P1.bound <- object$P1.bound
      }

      if (!is.null(object$P2.bound)) {
        out$measures$P2.bound <- object$P2.bound
      }
    }
  }
  out
}


#' @param x an object of class "summary.excurobj", usually, a result of a call
#'   to [summary.excurobj()].
#' @export
#' @method print summary.excurobj
#' @rdname summary.excurobj
print.summary.excurobj <- function(x, ...) {
  cat("Call: \n")
  print(x$call)
  cat("\nComputation:\n")
  cat(x$computation, "\n\n")

  if (x$calculation == "excursions") {
    cat("Level: u = ")
    cat(x$u, "\n")
    cat("Error probability: alpha = ")
    cat(x$alpha, "\n")
    cat("Limit for excursion function computation: F.limit = ")
    cat(x$F.limit, "\n")
    cat("Method used : ")
    cat(x$method, "\n")
  } else if (x$calculation == "simconf") {
    cat("Error probability: alpha = ")
    cat(x$alpha, "\n")
  } else if (x$calculation == "contourmap") {
    cat("Level: u = ")
    cat(x$u, "\n")
    cat("Type of contour map: ")
    cat(x$type, "\n")
    if (x$F.computed) {
      cat("Contour map function computed\n")
      cat("Limit for excursion function computation: F.limit = ")
      cat(x$F.limit, "\n")
    } else {
      cat("Contour map function not computed\n")
    }
    cat("Quality measures computed : ")
    if (is.null(x$measures)) {
      cat("none\n")
    } else {
      for (i in seq_along(x$measures)) {
        cat(names(x$measures)[i], " = ", x$measures[[i]], "\n")
      }
    }
  }
}

#' @export
#' @method print excurobj
#' @rdname summary.excurobj
print.excurobj <- function(x, ...) {
  print(summary(x))
}



# .onLoad <- function(libname, pkgname) {
#  # For Matrix coercion deprecation testing: 1=warn, 2=stop, NA=something else
#  # options(Matrix.warnDeprecatedCoerce = 2)
# }
