## For a result object computed in experimental mode, add back the linear
## predictor to the configs
inla.add.linearpredictor <- function(result, ind = NULL) {
  tau <- 1e9
  A <- rbind(
    result$misc$configs$pA %*% result$misc$configs$A,
    result$misc$configs$A
  )
  if (!is.null(ind)) {
    A <- A[ind, ]
  }

  I <- Diagonal(dim(A)[1])
  Abar <- rbind(cbind(I, -A), cbind(-t(A), t(A) %*% A))
  for (i in 1:result$misc$configs$nconfig) {
    Q <- bdiag(
      Matrix(0, nrow = dim(A)[1], ncol = dim(A)[1]),
      result$misc$configs$config[[i]]$Q
    )
    result$misc$configs$config[[i]]$Q <- Q + tau * Abar

    result$misc$configs$config[[i]]$Qinv <- bdiag(
      A %*% result$misc$configs$config[[i]]$Qinv %*% t(A),
      result$misc$configs$config[[i]]$Qinv
    )

    result$misc$configs$config[[i]]$mean <- c(
      as.double(A %*% result$misc$configs$config[[i]]$mean),
      result$misc$configs$config[[i]]$mean
    )
    result$misc$configs$config[[i]]$improved.mean <- c(
      as.double(A %*% result$misc$configs$config[[i]]$improved.mean),
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
        index <- 1:length(index.original)
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
        index <- 1:length(index.original)
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

private.get.config <- function(result, i) {
  mu <- result$misc$configs$config[[i]]$mean
  Q <- forceSymmetric(result$misc$configs$config[[i]]$Q)
  vars <- diag(result$misc$configs$config[[i]]$Qinv)
  m <- max(unlist(lapply(
    result$misc$configs$config,
    function(x) x$log.posterior
  )))
  lp <- result$misc$configs$config[[i]]$log.posterior - m
  return(list(mu = mu, Q = Q, vars = vars, lp = lp))
}

## Calculate the marginal probability for X_i>u or X_i<u.
## Note that the index 'i' refers to a location in the linear
## predictor if predictor==TRUE, whereas it refers to a location
## in the random effect vector otherwise.
inla.get.marginal <- function(i, u, result, effect.name = NULL, u.link, type) {
  if (is.null(effect.name) && u.link == TRUE) {
    marg.p <- result$marginals.fitted.values[[i]]
  } else if (is.null(effect.name)) {
    # Calculate marginals using linear predictor
    marg.p <- result$marginals.linear.predictor[[i]]
  } else {
    # Calculate marginals using a random effect
    marg.p <- result$marginals.random[[effect.name]][[i]]
  }

  if (type == "<") {
    return(INLA::inla.pmarginal(u, marg.p))
  } else {
    return(1 - INLA::inla.pmarginal(u, marg.p))
  }
}

## Calculate the marginal probability for a<X_i<b.
## The function returns c(P(X<a),P(X<b))
inla.get.marginal.int <- function(i, a, b, result, effect.name = NULL) {
  if (is.null(effect.name)) {
    # Calculate marginals using linear predictor
    marg.p <- result$marginals.linear.predictor[[i]]
  } else {
    # Calculate marginals using a random effect
    marg.p <- result$marginals.random[[effect.name]][[i]]
  }
  return(c(INLA::inla.pmarginal(a, marg.p), INLA::inla.pmarginal(b, marg.p)))
}
