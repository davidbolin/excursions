## Functions that will soon appear in the INLA package.

#library(INLA)

inla.mesh.segment <- function(...) {
    UseMethod("inla.mesh.segment")
}

inla.mesh.segment.default <-
    function(loc = NULL, idx = NULL, grp = NULL, is.bnd = TRUE, ...)
{
    if ((missing(loc) || is.null(loc)) &&
        (missing(idx) || is.null(idx))) {
        stop("At most one of 'loc' and 'idx' may be missing or null.")
    }
    if (!missing(loc) && !is.null(loc)) {
        if (!is.matrix(loc)) {
            loc = as.matrix(loc)
        }
        if (!is.double(loc)) {
            storage.mode(loc) = "double"
        }
        if (missing(idx) || is.null(idx))
            idx = (inla.ifelse(is.bnd,
                               c(1:nrow(loc),1),
                               c(1:nrow(loc))))
    } else {
        loc = NULL
    }

    if (!missing(idx) && !is.null(idx)) {
        if (!is.vector(idx) && !is.matrix(idx))
            stop("'idx' must be a vector or a matrix")
        if (is.vector(idx))
            idx = as.matrix(idx, nrow=length(idx), ncol=1)
        if (ncol(idx) == 1) {
            if (nrow(idx) < 2) {
                if (nrow(idx) == 1) {
                    warning("Segment specification must have at least 2, or 0, indices.")
                }
                idx <- matrix(0L, 0, 2)
            } else {
                idx = matrix(c(idx[-nrow(idx)],idx[-1]), nrow(idx)-1, 2)
            }
        }
        storage.mode(idx) <- "integer"
        if (!is.null(loc) &&
            (nrow(idx) > 0) &&
            (max(idx, na.rm=TRUE) > nrow(loc))) {
            warning("Segment indices (max=", max(idx, na.rm=TRUE),
                    ") exceed specified location list length (",
                    nrow(loc), ").")
        }
    }

    if (!missing(grp) && !is.null(grp)) {
        if (!is.vector(grp) && !is.matrix(grp))
            stop("'grp' must be a vector or a matrix")
        grp = matrix(grp, min(length(grp), nrow(idx)), 1)
        if (nrow(grp)<nrow(idx))
            grp = (matrix(c(as.vector(grp),
                            rep(grp[nrow(grp)], nrow(idx)-length(grp))),
                            nrow(idx), 1))
        storage.mode(grp) <- "integer"
    } else
        grp = NULL

    ## Filter away NAs in loc and idx
    if (!is.null(loc)) {
        idx[is.na(idx)] = 0L ## Avoid R annoyances with logical+NA indexing
        while (sum(is.na(loc))>0) {
            i = min(which(rowSums(is.na(loc))>0))
            loc = loc[-i,,drop=FALSE]
            idx[idx==i] = 0L
            idx[idx>i] = idx[idx>i]-1L
        }
        idx[idx==0L] = NA
    }
    while (sum(is.na(idx))>0) {
        i = min(which(rowSums(is.na(idx))>0))
        idx = idx[-i,,drop=FALSE]
        if (!is.null(grp))
            grp = grp[-i,,drop=FALSE]
    }

    if (!is.null(loc)) {
        ## Identify unused locations and remap indices accordingly.
        idx.new = rep(0L, nrow(loc))
        idx.new[as.vector(idx)] = 1L
        loc = loc[idx.new==1L,, drop=FALSE]
        idx.new[idx.new==1L] = seq_len(sum(idx.new))
        idx = (matrix(idx.new[as.vector(idx)],
                      nrow=nrow(idx),
                      ncol=ncol(idx)))
    }

    ret = list(loc=loc, idx=idx, grp=grp, is.bnd=is.bnd)
    class(ret) <- "inla.mesh.segment"
    return(ret)
}

inla.mesh.segment.inla.mesh.segment <- function(..., grp.default=0) {
    segm <- list(...)
    if (!all(unlist(lapply(segm,
                           function(x) inherits(x,"inla.mesh.segment"))))) {
        stop("All objects must be of class 'inla.mesh.segment'.")
    }

    Nloc <- unlist(lapply(segm, function(x) nrow(x$loc)))
    cumNloc <- c(0, cumsum(Nloc))
    Nidx <- unlist(lapply(segm, function(x) nrow(x$idx)))

    loc <- do.call(rbind, lapply(segm, function(x) x$loc))
    idx <- do.call(rbind, lapply(seq_along(segm),
                                function(x) segm[[x]]$idx+cumNloc[x]))
    grp <- unlist(lapply(seq_along(segm),
                         function(x) {
                             if (is.null(segm[[x]]$grp)) {
                                 rep(grp.default, Nidx[x])
                             } else {
                                 segm[[x]]$grp
                             }
                         }))
    is.bnd <- unlist(lapply(segm, function(x) x$is.bnd))
    if (!all(is.bnd) || all(!is.bnd)) {
        warning("Inconsistent 'is.bnd' attributes.  Setting 'is.bnd=FALSE'.")
        is.bnd <- FALSE
    } else {
        is.bnd <- all(is.bnd)
    }

    inla.mesh.segment(loc=loc, idx=idx, grp=grp, is.bnd=is.bnd)
}
