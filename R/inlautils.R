## For a result object computed in experimental mode, add back the linear
## predictor to the configs
inla.add.linearpredictor <- function(result, ind = NULL) {
  ## Precision of the link between the linear predictor and the latent field.
  ## Uses the INLA default predictor precision; much larger values make the
  ## joint precision numerically singular when A has dense columns (e.g. an
  ## intercept), causing the Cholesky factorisation to fail.
  tau <- exp(15)
  A <- rbind(
    result$misc$configs$pA %*% result$misc$configs$A,
    result$misc$configs$A
  )
  ## Offsets are not part of the latent field, so they must be added to the
  ## predictor mean. They are stored for the (A)Predictor rows only.
  offsets <- rep(0, dim(A)[1])
  off <- result$misc$configs$offsets
  if (length(off) > 0) {
    offsets[seq_along(off)] <- off
  }
  if (!is.null(ind)) {
    A <- A[ind, , drop = FALSE]
    offsets <- offsets[ind]
  }

  I <- Diagonal(dim(A)[1])
  Abar <- rbind(cbind(I, -A), cbind(-t(A), t(A) %*% A))
  for (i in 1:result$misc$configs$nconfig) {
    Q <- bdiag(
      Matrix(0, nrow = dim(A)[1], ncol = dim(A)[1]),
      result$misc$configs$config[[i]]$Q
    )
    result$misc$configs$config[[i]]$Q <- Q + tau * Abar

    ## The Qinv stored by INLA is only computed on the sparsity pattern of Q,
    ## so A %*% Qinv %*% t(A) misses covariances needed for the predictor
    ## variances. They are computed from the joint precision instead, so that
    ## they are consistent with the Q used in the integration. This is done in
    ## private.get.config(), only for the configurations that are used.
    result$misc$configs$config[[i]]$n.predictor <- dim(A)[1]

    result$misc$configs$config[[i]]$mean <- c(
      offsets + as.double(A %*% result$misc$configs$config[[i]]$mean),
      result$misc$configs$config[[i]]$mean
    )
    result$misc$configs$config[[i]]$improved.mean <- c(
      offsets + as.double(A %*% result$misc$configs$config[[i]]$improved.mean),
      result$misc$configs$config[[i]]$improved.mean
    )
  }
  result
}

## Find the indices into inla output config structures corresponding
## to a specific predictor, effect name, or inla.stack tag.
##
## result : an inla object
inla.output.indices <- function(result, name = NULL, stack = NULL, tag = NULL,
                                compressed = TRUE) {
  if (!is.null(result$misc$configs$.preopt) && result$misc$configs$.preopt) {
    inla.experimental <- TRUE
  } else {
    inla.experimental <- FALSE
  }

  if (!is.null(name) && !is.null(tag)) {
    stop("At most one of 'name' and 'tag' may be non-null.")
  }
  if (!is.null(tag)) {
    if (is.null(stack)) {
      stop("'tag' specified but 'stack' is NULL.")
    }
    tags <- names(stack$data$index)
    if (!(tag %in% tags)) {
      stop("'tag' not found in 'stack'.")
    }
  } else if (is.null(name) || (!is.null(name) && (name == ""))) {
    if ("APredictor" %in% result$misc$configs$contents$tag) {
      name <- "APredictor"
    } else {
      name <- "Predictor"
    }
  }
  result.updated <- FALSE

  ## Find variables
  if (!is.null(name)) {
    if (!(name %in% result$misc$configs$contents$tag)) {
      stop("'name' not found in result.")
    }
    ct <- result$misc$configs$contents

    # Shift indices for experimental mode
    if (inla.experimental && !(name %in% c("APredictor", "Predictor"))) {
      for (nm in c("APredictor", "Predictor")) {
        if (ct$tag[1] == nm) {
          ct$tag <- ct$tag[-1]
          ct$start <- ct$start[-1] - ct$start[2] + 1
          ct$length <- ct$length[-1]
        }
      }
      nameindex <- which(ct$tag == name)
      index <- (ct$start[nameindex] - 1L + seq_len(ct$length[nameindex]))
    } else if (inla.experimental) {
      # only add the part to be predicted
      nameindex <- which(ct$tag == name)
      index.original <- (ct$start[nameindex] - 1L + seq_len(ct$length[nameindex]))
      if (compressed) {
        result <- inla.add.linearpredictor(result, index.original)
        index <- seq_along(index.original)
      } else {
        result <- inla.add.linearpredictor(result)
        index <- index.original
      }
      result.updated <- TRUE
    } else {
      nameindex <- which(ct$tag == name)
      index <- (ct$start[nameindex] - 1L + seq_len(ct$length[nameindex]))
    }
  } else { ## Have tag
    index.original <- stack$data$index[[tag]]
    if (inla.experimental) {
      # only add the part to be predicted
      if (compressed) {
        result <- inla.add.linearpredictor(result, index.original)
        index <- seq_along(index.original)
      } else {
        result <- inla.add.linearpredictor(result)
        index <- index.original
      }
      result.updated <- TRUE
    } else {
      index <- index.original
    }
  }
  if (result.updated) {
    return(list(
      index = index,
      index.original = index.original,
      result = result,
      result.updated = result.updated
    ))
  } else {
    return(list(index = index, result.updated = result.updated))
  }
}

private.simconf.link <- function(res, links, trans = TRUE) {
  if (trans) {
    n <- length(res$a)
    res$a.marginal <- sapply(1:n, function(i) {
      private.link.function(
        res$a.marginal[i], links[i],
        inv = TRUE
      )
    })
    res$b.marginal <- sapply(1:n, function(i) {
      private.link.function(
        res$b.marginal[i], links[i],
        inv = TRUE
      )
    })
    res$a <- sapply(1:n, function(i) {
      private.link.function(
        res$a[i], links[i],
        inv = TRUE
      )
    })
    res$b <- sapply(1:n, function(i) {
      private.link.function(
        res$b[i], links[i],
        inv = TRUE
      )
    })
  }
  return(res)
}

private.link.function <- function(x, link, inv = FALSE) {
  if (is.na(link)) {
    link <- "identity"
  }
  return(do.call(paste("inla.link.", link, sep = ""), list(x = x, inv = inv)))
}

## Index of the configuration with the largest log posterior (the first one,
## if there are ties)
private.mode.config <- function(result) {
  which.max(vapply(
    result$misc$configs$config,
    function(x) x$log.posterior,
    0.0
  ))
}

private.get.config <- function(result, i) {
  mu <- result$misc$configs$config[[i]]$mean
  Q <- forceSymmetric(result$misc$configs$config[[i]]$Q)
  Qinv <- result$misc$configs$config[[i]]$Qinv
  vars <- diag(Qinv)
  n.pred <- result$misc$configs$config[[i]]$n.predictor
  if (!is.null(n.pred)) {
    ## Linear predictor added by inla.add.linearpredictor(). The Qinv of INLA
    ## is only for the latent field, so compute the selected inverse of the
    ## joint precision, and keep its covariances on the pattern of Q.
    sel <- private.selected.inverse(Q)
    vars <- c(sel$vars[seq_len(n.pred)], vars)
    Qt <- as(Matrix::triu(Q), "TsparseMatrix")
    Qinv <- Matrix::sparseMatrix(
      i = Qt@i + 1L, j = Qt@j + 1L,
      x = sel$cov(Qt@i + 1L, Qt@j + 1L), dims = dim(Q)
    )
  }
  m <- max(unlist(lapply(
    result$misc$configs$config,
    function(x) x$log.posterior
  )))
  lp <- result$misc$configs$config[[i]]$log.posterior - m

  ## Qinv has the covariances on (one triangle of) the pattern of Q
  list(mu = mu, Q = Q, vars = vars, Qinv = Qinv, lp = lp)
}

## The marginals of the nodes i of the linear predictor (with predictor ==
## TRUE), the fitted values (with u.link), or the random effect effect.name.
inla.get.marginals <- function(i, result, effect.name = NULL, u.link = FALSE) {
  if (is.null(effect.name) && u.link) {
    result$marginals.fitted.values[i]
  } else if (is.null(effect.name)) {
    result$marginals.linear.predictor[i]
  } else {
    result$marginals.random[[effect.name]][i]
  }
}

## Calculate the marginal probabilities for X_i>u or X_i<u for the nodes i.
## Note that the index 'i' refers to a location in the linear
## predictor if predictor==TRUE, whereas it refers to a location
## in the random effect vector otherwise.
inla.get.marginal <- function(i, u, result, effect.name = NULL, u.link, type) {
  p <- private.pmarginals(inla.get.marginals(i, result, effect.name, u.link), u)
  if (type == "<") {
    p
  } else {
    1 - p
  }
}

## Calculate the marginal probabilities for a<X_i<b for the nodes i. Returns
## the matrix with columns P(X_i<a) and P(X_i<b).
inla.get.marginal.int <- function(i, a, b, result, effect.name = NULL) {
  m <- inla.get.marginals(i, result, effect.name)
  cbind(private.pmarginals(m, a), private.pmarginals(m, b))
}

## The distribution functions of the marginals at q, which is recycled to the
## number of marginals, where each marginal is a matrix or list with the
## values x and the densities y. This computes the same as
## INLA::inla.pmarginal(q[k], marginals[[k]]) for each k, which interpolates
## the log density with a spline and integrates it numerically, but for all
## marginals together, which is much faster. The log density is interpolated
## by cubic Hermite polynomials, with the derivatives from the differences to
## the neighbouring points, and these are integrated with Gauss-Legendre
## quadrature in each interval. As for inla.pmarginal, the points with
## negligible density are removed (see INLA:::inla.marginal.fix), the
## distribution is normalised on the range of x, and q is truncated to it.
##
## The marginals of INLA are usually matrices with the same number of rows,
## which are then stacked at once, and only the marginals with points to
## remove are handled one at a time.
private.pmarginals <- function(marginals, q) {
  K <- length(marginals)
  q <- rep_len(q, K)
  p <- numeric(K)
  if (K == 0) {
    return(p)
  }
  eps <- .Machine$double.eps * 1000
  rest <- seq_len(K)
  len <- lengths(marginals, use.names = FALSE)
  if (all(len == len[1]) && len[1] >= 4 && len[1] %% 2 == 0 &&
    all(vapply(marginals, function(m) is.matrix(m) && ncol(m) == 2L, TRUE))) {
    np <- len[1] / 2
    A <- array(unlist(marginals, use.names = FALSE), c(np, 2L, K))
    X <- t(A[, 1L, ])
    Y <- t(A[, 2L, ])
    if (K == 1) {
      X <- matrix(X, 1)
      Y <- matrix(Y, 1)
    }
    ## The marginals where no points have to be removed
    ok <- !anyNA(Y)
    ok <- if (ok) rowSums(!(Y > 0)) == 0 else !apply(is.na(Y) | !(Y > 0), 1, any)
    ymax <- numeric(K)
    ymax[ok] <- Y[ok, , drop = FALSE][cbind(seq_len(sum(ok)), max.col(Y[ok, , drop = FALSE], "first"))]
    ok[ok] <- rowSums(Y[ok, , drop = FALSE] / ymax[ok] <= eps) == 0
    if (any(ok)) {
      p[ok] <- private.pmarginals.matrix(
        X[ok, , drop = FALSE], log(Y[ok, , drop = FALSE]), q[ok]
      )
    }
    rest <- which(!ok)
  }
  if (length(rest) == 0) {
    return(p)
  }
  xy <- lapply(marginals[rest], function(m) {
    if (is.matrix(m)) {
      x <- m[, 1L]
      y <- m[, 2L]
    } else {
      x <- m[["x"]]
      y <- m[["y"]]
    }
    ok <- !is.na(y)
    x <- x[ok]
    y <- y[ok]
    ok <- y > 0 & y / max(y) > eps
    list(x = x[ok], y = y[ok])
  })
  len <- vapply(xy, function(m) length(m$x), 1L)
  for (np in unique(len)) {
    k <- which(len == np)
    if (np < 2) {
      ## A single point, as a point mass
      x1 <- vapply(xy[k], function(m) if (np == 1) m$x else NA_real_, 0)
      p[rest[k]] <- as.numeric(q[rest[k]] >= x1)
      next
    }
    X <- matrix(unlist(lapply(xy[k], `[[`, "x")), ncol = np, byrow = TRUE)
    Y <- log(matrix(unlist(lapply(xy[k], `[[`, "y")), ncol = np, byrow = TRUE))
    p[rest[k]] <- private.pmarginals.matrix(X, Y, q[rest[k]])
  }
  p
}

## private.pmarginals for marginals with the values in the rows of X, at least
## two, and the log densities in the rows of Y.
private.pmarginals.matrix <- function(X, Y, q) {
  nk <- nrow(X)
  np <- ncol(X)
  ## Gauss-Legendre nodes and weights on [0, 1]
  gs <- c(-0.9061798459386640, -0.5384693101056831, 0, 0.5384693101056831, 0.9061798459386640)
  gw <- c(0.2369268850561891, 0.4786286704993665, 0.5688888888888889, 0.4786286704993665, 0.2369268850561891)
  gs <- (gs + 1) / 2
  gw <- gw / 2
  h <- X[, -1, drop = FALSE] - X[, -np, drop = FALSE]
  d <- (Y[, -1, drop = FALSE] - Y[, -np, drop = FALSE]) / h
  ## Derivatives of the log density at the points
  D <- matrix(0, nk, np)
  if (np == 2) {
    D[, 1] <- D[, 2] <- d[, 1]
  } else {
    h1 <- h[, -(np - 1), drop = FALSE]
    h2 <- h[, -1, drop = FALSE]
    D[, 2:(np - 1)] <- (h2 * d[, -(np - 1), drop = FALSE] + h1 * d[, -1, drop = FALSE]) / (h1 + h2)
    D[, 1] <- ((2 * h[, 1] + h[, 2]) * d[, 1] - h[, 1] * d[, 2]) / (h[, 1] + h[, 2])
    D[, np] <- ((2 * h[, np - 1] + h[, np - 2]) * d[, np - 1] - h[, np - 1] * d[, np - 2]) /
      (h[, np - 1] + h[, np - 2])
  }
  ## The log density on [0, 1] of each interval is the cubic Hermite
  ## polynomial with the values y0, y1 and the scaled derivatives m0, m1.
  ## Subtract the maximum before exponentiating, for the tails.
  ymax <- Y[cbind(seq_len(nk), max.col(Y, "first"))]
  y0 <- Y[, -np, drop = FALSE] - ymax
  y1 <- Y[, -1, drop = FALSE] - ymax
  m0 <- D[, -np, drop = FALSE] * h
  m1 <- D[, -1, drop = FALSE] * h
  ## Integral over [0, tau] of each interval, on the scale of x
  seg <- function(tau, y0, y1, m0, m1, h) {
    out <- 0
    for (g in seq_along(gs)) {
      t <- tau * gs[g]
      t2 <- t * t
      t3 <- t2 * t
      lf <- (2 * t3 - 3 * t2 + 1) * y0 + (t3 - 2 * t2 + t) * m0 +
        (-2 * t3 + 3 * t2) * y1 + (t3 - t2) * m1
      out <- out + gw[g] * exp(lf)
    }
    out * tau * h
  }
  full <- seg(1, y0, y1, m0, m1, h)
  total <- rowSums(full)
  ## The interval j of q, which is truncated to the range of x, the integral
  ## of the intervals before it, and of the part of it before q
  qk <- pmin(pmax(q, X[, 1]), X[, np])
  j <- pmin(rowSums(X <= qk), np - 1)
  before <- rowSums(full * (col(full) < j))
  ij <- cbind(seq_len(nk), j)
  tau <- (qk - X[ij]) / h[ij]
  part <- seg(tau, y0[ij], y1[ij], m0[ij], m1[ij], h[ij])
  pmin(pmax((before + part) / total, 0), 1)
}
