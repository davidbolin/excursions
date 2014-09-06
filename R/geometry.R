library(sp)
library(rgeos)
library(INLA)


## Trace a length 2 contour segment through pixels
##
contour.segment.pixels <- function(c.x, c.y,
                                   edge.x, edge.y,
                                   pairs=FALSE)
{
    nx <- length(edge.x)
    ny <- length(edge.y)
    d.x <- diff(c.x)
    d.y <- diff(c.y)
    if (abs(d.x) < abs(d.y)) {
        ## Swap x and y
        idx <- contour.segment.pixels(c.y, c.x, edge.y, edge.x, TRUE)
        idx <- list(x=idx$y, y=idx$x)
    } else if (d.y < 0) {
        ## Flip ordering of segment
        idx <- contour.segment.pixels(c.x[2:1], c.y[2:1], edge.x, edge.y, TRUE)
        idx <- list(x=idx$x[length(idx$x):1], y=idx$y[length(idx$y):1])
    } else {
        idx <- list(x=numeric(0), y=numeric(0))
        ## |d.x| >= d.y >= 0
        dir.x <- sign(d.x) ## Can be 0 only if dir.y is 0
        dir.y <- sign(d.y) ## Must now be >= 0
        if (dir.y == 0) {
            start.y <- min(which((c.y[1] >= edge.y[-ny]) &
                                 (c.y[1] <= edge.y[-1])))
            if (dir.x == 0) {
                start.x <- min(which((c.x[1] >= edge.x[-nx]) &
                                     (c.x[1] <= edge.x[-1])))
                end.x <- start.x
            } else if (dir.x > 0) {
                start.x <- min(which((c.x[1] >= edge.x[-nx]) &
                                     (c.x[1] < edge.x[-1])))
                end.x <- min(which((c.x[2] > edge.x[-nx]) &
                                   (c.x[2] <= edge.x[-1])))
            } else { ## dir.x < 0
                start.x <- min(which((c.x[1] > edge.x[-nx]) &
                                     (c.x[1] <= edge.x[-1])))
                end.x <- min(which((c.x[2] >= edge.x[-nx]) &
                                   (c.x[2] < edge.x[-1])))
            }
            idx$x <- start.x:end.x
            idx$y <- rep(start.y, length(idx$x))
        } else { ## dir.y > 0
            start.y <- min(which((c.y[1] >= edge.y[-ny]) &
                                 (c.y[1] < edge.y[-1])))
            end.y <- min(which((c.y[2] > edge.y[-ny]) &
                               (c.y[2] <= edge.y[-1])))
            if (dir.x > 0) {
                start.x <- min(which((c.x[1] >= edge.x[-nx]) &
                                     (c.x[1] < edge.x[-1])))
                end.x <- min(which((c.x[2] > edge.x[-nx]) &
                                   (c.x[2] <= edge.x[-1])))
            } else { ## dir.x < 0
                start.x <- min(which((c.x[1] > edge.x[-nx]) &
                                     (c.x[1] <= edge.x[-1])))
                end.x <- min(which((c.x[2] >= edge.x[-nx]) &
                                   (c.x[2] < edge.x[-1])))
            }
            if (start.y == end.y) {
                idx$x <- start.x:end.x
                idx$y <- rep(start.y, length(idx$x))
            } else {
                ## Need to step thorugh each y-pixel-row.
                curr.y <- start.y
                curr.x <- start.x
                while (curr.y < end.y) {
                    next.y <- curr.y + 1
                    ## Find intersection with edge.y[next.y]
                    ## (x-c.x[1])*d.y/d.x + c.y[1] = edge.y[next.y]
                    intersect.x <- c.x[1]+d.x/d.y*(edge.y[next.y]-c.y[1])
                    if (dir.x > 0) {
                        next.x <- min(which((intersect.x >= edge.x[-nx]) &
                                            (intersect.x < edge.x[-1])))
                    } else { ## dir.x < 0
                        next.x <- min(which((intersect.x > edge.x[-nx]) &
                                            (intersect.x <= edge.x[-1])))
                    }
                    idx$x <- c(idx$x, curr.x:next.x)
                    idx$y <- c(idx$y, rep(curr.y, length(curr.x:next.x)))
                    curr.y <- next.y
                    curr.x <- next.x
                }
                idx$x <- c(idx$x, curr.x:end.x)
                idx$y <- c(idx$y, rep(curr.y, length(curr.x:end.x)))
            }
        }
    }
    if (pairs) {
        return(idx)
    } else {
        return((idx$y - 1) * (nx - 1) + idx$x)
    }
}


##' Trace line segments through grid pixels
##'
##' .. content for \details{} ..
##' @title Discretise line segments
##' @param contourlines contourLines
##' @param pixelgrid 
##' @param do.plot 
##' @return Pixel index vector
##' @author Finn Lindgren
contour.pixels <- function(contourlines, pixelgrid, do.plot=0)
{
    if (inherits(pixelgrid, "pgrid")) {
        mid.x <- pixelgrid$upx
        mid.y <- pixelgrid$upy
        edge.x <- pixelgrid$ubx
        edge.y <- pixelgrid$uby
    } else if (inherits(pixelgrid, "inla.mesh.lattice")) {
        mid.x <- pixelgrid$x
        mid.y <- pixelgrid$y
        step.x <- mid.x[2]-mid.x[1]
        step.y <- mid.y[2]-mid.y[1]
        edge.x <- c(mid.x[1]-step.x/2, mid.x+step.x/2)
        edge.y <- c(mid.y[1]-step.y/2, mid.y+step.y/2)
    } else {
        stop("Unsupported grid specification class.")
    }
    nx <- length(mid.x)
    ny <- length(mid.y)
    which.pixel <- vector("list", length(contourlines))
    for (level in seq_along(contourlines)) {
        c.x <- contourlines[[level]]$x
        c.y <- contourlines[[level]]$y
        for (segment in seq_len(length(c.x)-1)) {
            tmp <- contour.segment.pixels(c.x[segment+(0:1)],
                                          c.y[segment+(0:1)],
                                          edge.x,
                                          edge.y)
            which.pixel[[level]] <- c(which.pixel[[level]], tmp)
            if (do.plot > 1) {
                setvec <- rep(0, nx*ny)
                setvec[sort(unique(unlist(which.pixel)))] <- 1
                image(mid.x, mid.y, matrix(setvec, nx, ny), col = c(0, 3))
                setvec <- rep(0, nx*ny)
                setvec[tmp] <- 1
                image(mid.x, mid.y, matrix(setvec, nx, ny), col = c(0, 4),
                      add=TRUE)
                plot.contourLines(contourlines, add=TRUE)
                lines(c.x[segment+(0:1)], c.y[segment+(0:1)], col=2)
                if (do.plot > 2) {
                    readline("next")
                }
            }
        }
    }
    if (do.plot > 0) {
        setvec <- rep(0, nx*ny)
        setvec[sort(unique(unlist(which.pixel)))] <- 1
        image(mid.x, mid.y, matrix(setvec, nx, ny), col = c(0, 3))
        plot.contourLines(contourlines, add=TRUE)
        lines(c.x[segment+(0:1)], c.y[segment+(0:1)], col=2)
    }

    sort(unique(unlist(which.pixel)))
}


## Connect segments into sequences looping around regions,
## in counterclockwise order.
## Initial segments have inner on their left hand side.
## grp.ccw keeps its orientation (multiple groups allowed)
## grp.cw is traversed in reverse orientation (multiple groups allowed)
##
## while (any segments remain) do
##   while (forward connected segments are found) do
##     follow segments forward
##   if (sequence not closed loop)
##     while (backward connected segments are found) do
##       follow segments backward
##
## sequences : list
## sequences[[k]] : node index vector for a single sequence
## seg : list
## seg[[k]] : segment index vector for a single sequence
## grp : list
## grp[[k]] : group index vector for a single sequence
connect.segments <-function(segment.set,
                            segment.grp=rep(0L, nrow(segment.set)),
                            grp.ccw=unique(segment.grp),
                            grp.cw=integer(0),
                            ccw=TRUE,
                            ambiguous.warning=FALSE)
{
    ## Remove unneeded segments
    segment.idx <- seq_len(nrow(segment.set))
    segment.idx <- c(segment.idx[segment.grp %in% grp.ccw],
                     segment.idx[segment.grp %in% grp.cw])
    segment.set <- segment.set[segment.idx,,drop=FALSE]
    segment.grp <- as.integer(segment.grp[segment.idx])
    ## Unneeded segment removal done
    ## Reverse the direction of segments in grp.cw
    segment.set[segment.grp %in% grp.cw, ] <-
        segment.set[segment.grp %in% grp.cw, 2:1, drop=FALSE]
    ## Segment reversion done
    nE <- nrow(segment.set)
    if (nE == 0) {
        return(list(sequences=list(), seg=list(), grp=list()))
    }
    ## Remap nodes into 1...nV
    segments <- sort(unique(as.vector(segment.set)))
    nV <- length(segments)
    segments.reo <- spam(list(i=segments,
                              j=rep(1, nV),
                              values=seq_len(nV)),
                         max(segments), 1)
    segment.set <- matrix(segments.reo[as.vector(segment.set)],
                          nrow(segment.set), ncol(segment.set))
    ## Node remapping done

    segment.unused <- rep(TRUE, nE)
    segment.unused.idx <- which(segment.unused)
    segment.VV <- spam(list(i=segment.set[,1],
                            j=segment.set[,2],
                            values=seq_len(nE)),
                       nV, nV)
    loops.seg <- list()
    loops <- list()
    grp <- list()
    while (length(segment.unused.idx) > 0) {
        n <- 0
        loop.seg <- integer(0)
        ## Forward loop
##        message("Forwards")
        while (length(segment.unused.idx) > 0) {
            if (n == 0) {
                si <- segment.unused.idx[1]
            } else {
                si <- as.vector(as.matrix(segment.VV[segment.set[si,2],]))
                si <- si[si %in% segment.unused.idx]
                if (length(si) == 0) {
                    ## End of sequence
                    break
                } else {
                    if ((length(si) > 1) && ambiguous.warning) {
                        warning("Ambiguous segment sequence.")
                    }
                    si <- si[1]
                }
            }
            segment.unused.idx <-
                segment.unused.idx[segment.unused.idx != si]
            segment.unused[si] <- FALSE
            loop.seg <- c(loop.seg, si)
            n <- n+1

##            print(loop.seg)
##            loop <- c(segment.set[loop.seg,1],
##                      segment.set[loop.seg[n],2])
##            print(loop)
        }
        if ((segment.set[loop.seg[n],2] != segment.set[loop.seg[1],1]) &&
            (length(segment.unused.idx) > 0)) {
            ## No closed sequence found
            ## Backward loop
##            message("Backwards")
            si <- loop.seg[1]
            while (length(segment.unused.idx) > 0) {
                si <- as.vector(as.matrix(segment.VV[,segment.set[si,1]]))
                si <- si[si %in% segment.unused.idx]
                if (length(si) == 0) {
                    ## End of sequence
                    break
                } else if (length(si) > 1) {
                    warning("Ambiguous segment sequence.")
                    si <- si[1]
                } else {
                    si <- si[1]
                }
                segment.unused.idx <-
                    segment.unused.idx[segment.unused.idx != si]
                segment.unused[si] <- FALSE
                loop.seg <- c(si, loop.seg)
                n <- n+1

##                print(loop.seg)
##                loop <- c(segment.set[loop.seg,1],
##                          segment.set[loop.seg[n],2])
##                print(loop)
            }
        }

        loop <- c(segment.set[loop.seg,1],
                  segment.set[loop.seg[n],2])
        loop.grp <- segment.grp[loop.seg]

        loops.seg <- c(loops.seg, list(loop.seg))
        loops <- c(loops, list(loop))
        grp <- c(grp, list(loop.grp))
    }

    ## Remap nodes and segments back to original indices
    for (k in seq_along(loops)) {
        loops[[k]] <- segments[loops[[k]]]
        loops.seg[[k]] <- segment.idx[loops.seg[[k]]]
    }
    ## Node and segment index remapping done

    if (!ccw) {
        for (k in seq_along(loops)) {
            loops[[k]] <- loops[[k]][length(loops[[k]]):1]
            loops.seg[[k]] <- loops.seg[[k]][length(loops.seg[[k]]):1]
            grp[[k]] <- grp[[k]][length(grp[[k]]):1]
        }
    }

    return(list(sequences=loops, seg=loops.seg, grp=grp))
}

## Compute simple outline of 1/0 set on a grid, eliminating spikes.
outline.on.grid <- function(z, x=NULL, y=NULL)
{
    if (!require(spam)) {
        stop("The 'spam' package is needed.")
    }
    ni <- nrow(z)
    nj <- ncol(z)
    z <- (z != FALSE)
    if (missing(x) || is.null(x)) {
        x <- seq(0,1,length=ni)
    }
    if (missing(y) || is.null(y)) {
        y <- seq(0,1,length=nj)
    }

    ij2k <-function(i,j) {
        return((j-1)*ni+i)
    }

    seg <- matrix(integer(), 0, 2)
    bnd.seg <- matrix(integer(), 0, 2)

    ## Extract horizontal segment locations:
    zz.p <- z[-ni,-c(1,2),drop=FALSE] + z[-1,-c(1,2),drop=FALSE]
    zz.n <- z[-ni,-c(1,nj),drop=FALSE] + z[-1,-c(1,nj),drop=FALSE]
    zz.m <- z[-ni,-c(nj-1,nj),drop=FALSE] + z[-1,-c(nj-1,nj),drop=FALSE]
    zz <- (zz.n == 2) * (zz.p+zz.m > 0) * ((zz.p == 0)*1 + (zz.m == 0)*2)
    ## zz=0 : No segment, zz=1 : set is below, zz=2 : set is above
    ijv <- triplet(as.spam(zz == 1), tri=TRUE)
    idx <- which(ijv$values > 0)
    seg <- rbind(seg,
                 cbind(ij2k(ijv$i[idx]+1, ijv$j[idx]+1),
                       ij2k(ijv$i[idx], ijv$j[idx]+1)))
    ijv <- triplet(as.spam(zz == 2), tri=TRUE)
    idx <- which(ijv$values > 0)
    seg <- rbind(seg,
                 cbind(ij2k(ijv$i[idx], ijv$j[idx]+1),
                       ij2k(ijv$i[idx]+1, ijv$j[idx]+1)))

    ## Extract vertical segment locations:
    zz.p <- z[-c(1,2),-nj,drop=FALSE] + z[-c(1,2),-1,drop=FALSE]
    zz.n <- z[-c(1,ni),-nj,drop=FALSE] + z[-c(1,ni),-1,drop=FALSE]
    zz.m <- z[-c(ni-1,ni),-nj,drop=FALSE] + z[-c(ni-1,ni),-1,drop=FALSE]
    zz <- (zz.n == 2) * (zz.p+zz.m > 0) * ((zz.p == 0)*1 + (zz.m == 0)*2)
    ## zz=0 : No segment, zz=1 : set is on left, zz=2 : set is on right
    ijv <- triplet(as.spam(zz == 1), tri=TRUE)
    idx <- which(ijv$values > 0)
    seg <- rbind(seg,
                 cbind(ij2k(ijv$i[idx]+1, ijv$j[idx]),
                       ij2k(ijv$i[idx]+1, ijv$j[idx]+1)))
    ijv <- triplet(as.spam(zz == 2), tri=TRUE)
    idx <- which(ijv$values > 0)
    seg <- rbind(seg,
                 cbind(ij2k(ijv$i[idx]+1, ijv$j[idx]+1),
                       ij2k(ijv$i[idx]+1, ijv$j[idx])))

    ## Extract diagonal segment locations;
    ## Quadruples with three 1, one 0
    zz <- z[-ni,-nj,drop=FALSE] + z[-1,-nj,drop=FALSE] + z[-ni,-1,drop=FALSE] + z[-1,-1,drop=FALSE]
    ## Which element was 0?
    zz <- (zz == 3) *
        (15 - (z[-ni,-nj,drop=FALSE]*1 + z[-1,-nj,drop=FALSE]*2 + z[-ni,-1,drop=FALSE]*4 + z[-1,-1,drop=FALSE]*8))
    ## zz=0 : No diagonal
    ## zz=1 : (0,0), zz=2 : (1,0), zz=4 : (0,1), zz=8 : (1,1)
    ijv <- triplet(as.spam(zz == 1), tri=TRUE)
    idx <- which(ijv$values > 0)
    seg <- rbind(seg,
                 cbind(ij2k(ijv$i[idx], ijv$j[idx]+1),
                       ij2k(ijv$i[idx]+1, ijv$j[idx])))
    ijv <- triplet(as.spam(zz == 2), tri=TRUE)
    idx <- which(ijv$values > 0)
    seg <- rbind(seg,
                 cbind(ij2k(ijv$i[idx], ijv$j[idx]),
                       ij2k(ijv$i[idx]+1, ijv$j[idx]+1)))
    ijv <- triplet(as.spam(zz == 4), tri=TRUE)
    idx <- which(ijv$values > 0)
    seg <- rbind(seg,
                 cbind(ij2k(ijv$i[idx]+1, ijv$j[idx]+1),
                       ij2k(ijv$i[idx], ijv$j[idx])))
    ijv <- triplet(as.spam(zz == 8), tri=TRUE)
    idx <- which(ijv$values > 0)
    seg <- rbind(seg,
                 cbind(ij2k(ijv$i[idx]+1, ijv$j[idx]),
                       ij2k(ijv$i[idx], ijv$j[idx]+1)))

    ## Extract horizontal boundary segment locations:
    zz.pm <- z[-ni,c(2,nj-1),drop=FALSE] + z[-1,c(2,nj-1),drop=FALSE]
    zz.n <- z[-ni,c(1,nj),drop=FALSE] + z[-1,c(1,nj),drop=FALSE]
    zz <- (zz.n == 2) * (zz.pm > 0) * matrix(rep(c(2,1), each=ni-1), ni-1, 2)
    ## zz=0 : No segment, zz=1 : set is below, zz=2 : set is above
    ijv <- triplet(as.spam(zz == 1), tri=TRUE)
    idx <- which(ijv$values > 0)
    bnd.seg <- rbind(bnd.seg,
                     cbind(ij2k(ijv$i[idx]+1, nj),
                           ij2k(ijv$i[idx], nj)))
    ijv <- triplet(as.spam(zz == 2), tri=TRUE)
    idx <- which(ijv$values > 0)
    bnd.seg <- rbind(bnd.seg,
                     cbind(ij2k(ijv$i[idx], 1),
                           ij2k(ijv$i[idx]+1, 1)))

    ## Extract vertical boundary segment locations:
    zz.pm <- z[c(2,ni-1),-nj,drop=FALSE] + z[c(2,ni-1),-1,drop=FALSE]
    zz.n <- z[c(1,ni),-nj,drop=FALSE] + z[c(1,ni),-1,drop=FALSE]
    zz <- (zz.n == 2) * (zz.pm > 0) * matrix(rep(c(2,1), times=nj-1), 2, nj-1)
    ## zz=0 : No segment, zz=1 : set is on left, zz=2 : set is on right
    ijv <- triplet(as.spam(zz == 1), tri=TRUE)
    idx <- which(ijv$values > 0)
    bnd.seg <- rbind(bnd.seg,
                     cbind(ij2k(ni, ijv$j[idx]),
                           ij2k(ni, ijv$j[idx]+1)))
    ijv <- triplet(as.spam(zz == 2), tri=TRUE)
    idx <- which(ijv$values > 0)
    bnd.seg <- rbind(bnd.seg,
                     cbind(ij2k(1, ijv$j[idx]+1),
                           ij2k(1, ijv$j[idx])))

    segment.grp <- rep(c(1L, 0L), c(nrow(seg), nrow(bnd.seg)))
    segment.set <- rbind(seg, bnd.seg)

    loc <- cbind(rep(x, times=nj), rep(y, each=ni))

    return(list(loc=loc, idx=segment.set, grp=segment.grp))
}

outline.to.sp <- function(outline,
                          grp.ccw=unique(outline$grp),
                          grp.cw=integer(0),
                          ccw=FALSE,
                          ambiguous.warning=FALSE,
                          ID="outline",
                          closed=TRUE,
                          ...)
{
    seg <- connect.segments(outline$idx, outline$grp,
                            grp.ccw=grp.ccw,
                            grp.cw=grp.cw,
                            ccw=ccw,
                            ambiguous.warning=ambiguous.warning)

    coords <- list()
    for (k in seq_along(seg$sequences)) {
        coords <- c(coords, list(outline$loc[seg$sequences[[k]],
                                             1:2, drop=FALSE]))
    }

    if (closed) {
        as.SpatialPolygons.raw(coords, ID=ID)
    } else {
        as.SpatialLines.raw(coords, ID=ID)
    }
}


outline.to.inla.mesh.segment <- function(outline,
                                         grp.ccw=unique(outline$grp),
                                         grp.cw=integer(0),
                                         ...)
{
   ik.ccw = outline$grp %in% grp.ccw
   ik.cw = outline$grp %in% grp.cw
   inla.mesh.segment(loc=outline$loc,
                     idx=rbind(outline$idx[ik.ccw,],
                               outline$idx[ik.cw, 2:1]),
                     grp=c(outline$grp[ik.ccw], outline$grp[ik.cw]))
}



as.SpatialPolygons.raw <- function(sequences, ID=" ") {
    polys <- lapply(sequences,
                    function(x) {
                        if (is.list(x)) {
                            p <- Polygon(cbind(x$x, x$y))
                        } else {
                            p <- Polygon(x)
                        }
                        p
                    })
    if (length(polys) == 0) {
        sp <- NULL
    } else {
        sp <- SpatialPolygons(list(Polygons(polys, ID=ID)))
    }
    sp
}

as.SpatialLines.raw <- function(cl, ID=" ") {
    polys <- lapply(cl,
                    function(x) {
                        if (is.list(x)) {
                            p <- Line(cbind(x$x, x$y))
                        } else {
                            p <- Line(x)
                        }
                        p
                    })
    if (length(polys) == 0) {
        sp <- NULL
    } else {
        sp <- SpatialLines(list(Lines(polys, ID=ID)))
    }
    sp
}





tricontour <- function(x, z, nlevels = 10,
                       levels = pretty(range(z, na.rm = TRUE), nlevels),
                       ...)
{
    message("tricontour")
    UseMethod("tricontour")
}

tricontour.inla.mesh <- function(x, z, nlevels = 10,
                       levels = pretty(range(z, na.rm = TRUE), nlevels),
                                 ...)
{
    message("tricontour.inla.mesh")
    tricontour.list(x$graph, z=z,
                    nlevels=nlevels, levels=levels,
                    loc=x$loc, ...)
}

tricontour.matrix <- function(x, z, nlevels = 10,
                              levels = pretty(range(z, na.rm = TRUE), nlevels),
                              loc, ...)
{
    message("tricontour.matrix")
    tricontour.list(list(tv=x), z=z,
                    nlevels=nlevels, levels=levels,
                    loc=loc, ...)
}



## Returns val=list(loc, idx, grp), where
##   grp = 1,...,nlevels*2+1, level groups are even, 2,4,...
## Suitable for
##   inla.mesh.segment(val$loc, val$idx[val$grp==k], val$idx[val$grp==k])
##     (supports R2 and S2)
## and, for odd k=1,3,...,nlevels*2-1,nlevels*2+1,
##   seg <- outline.to.inla.mesh.segment(val, grp.ccw=c(k-1,k), grp.cw=c(k+1))
##   sp <- outline.to.sp(val, grp.ccw=c(k-1,k), grp.cw=c(k+1), ccw=FALSE)
display.dim.list <- function(x) {
    lapply(as.list(sort(names(x))),
           function(xx) {
               sz <- dim(x[[xx]])
               type <- mode(x[[xx]])
               cl <- class(x[[xx]])[[1]]
               if (is.null(sz)) {
                   sz <- length(x[[xx]])
               }
               message(paste(xx, " = ", paste(sz, collapse=" x "),
                             " (", type, ", ", cl, ")", sep=""))
           }
           )
    invisible()
}
tricontour.list <- function(x, z, nlevels = 10,
                            levels = pretty(range(z, na.rm = TRUE), nlevels),
                            loc, type=c("+", "-"), tol=1e-7, ...)
{
    message("tricontour.list")
    type <- match.arg(type)
    nlevels <- length(levels)

    ## Generate graph properties
    x$Nt <- nrow(x$tv)
    x$Ne <- 3*x$Nt
    x$Nt <- max(as.vector(x$tv))
    x <- generate.graph.properties(x)
    x$ev <- cbind(as.vector(x$tv[,c(2,3,1)]),
                  as.vector(x$tv[,c(3,1,2)]))
    x$et <- rep(seq_len(x$Nt), times=3)
    x$eti <- rep(1:3, each=x$Nt) ## Opposing vertex within-triangle-indices
    x$te <- matrix(seq_len(x$Ne), x$Nt, 3)
    if (is.null(x$tt)) {
        stop("TODO: generate missing graph property 'tt'")
    }
    if (is.null(x$tti)) {
        stop("TODO: generate missing graph property 'tti'")
    }
    Nv <- x$Nv

    ## Find vertices on levels
    ## For each edge on a level, store edge if
    ##     opposing vertices are +/-, and either
    ##     0/- (type="+", u1 <= z < u2) or
    ##     +/0 (type="-", u1 < z <= u2)
    ## For each edge crossing at least one level,
    ##   calculate splitting vertices
    ##   if boundary edge, store new split edges
    ## For each triangle, find non-level edge crossings, and
    ##   store new vertex-edge crossing edges
    ##   store new edge-edge crossing edges
    idx <- matrix(0,0,2)
    grp <- integer(0)

    ## Find vertices on levels
    vcross.lev <- integer(length(z))
    for (lev in seq_along(levels)) {
        signv <- (z > levels[lev]+tol) - (z < levels[lev]-tol)
        vcross.lev[ signv == 0 ] <- lev
    }
    ## Find level crossing span for each edge (includes flat edges in levels)
    ecross.grp.lower <- rep(1L, x$Ne)
    ecross.grp.upper <- rep(2L*length(levels)+1L, x$Ne)
    for (lev in seq_along(levels)) {
        signv <- (z > levels[lev]+tol) - (z < levels[lev]-tol)
        lev.grp <- 2L*lev
        i <- pmin(signv[x$ev[,1]], signv[x$ev[,2]])
        ecross.grp.lower[ i == 0 ] <- lev.grp
        ecross.grp.lower[ i > 0 ] <- lev.grp+1L
    }
    for (lev in seq_along(levels)[length(levels):1L]) {
        signv <- (z > levels[lev]+tol) - (z < levels[lev]-tol)
        lev.grp <- 2L*lev
        i <- pmax(signv[x$ev[,1]], signv[x$ev[,2]])
        ecross.grp.upper[ i == 0 ] <- lev.grp
        ecross.grp.upper[ i  < 0 ] <- lev.grp-1L
    }

    ## For each edge on a level, store edge if ...
    ##   opposing vertices are +/-, and either
    ##   0/- (type="+", u1 <= z < u2) or
    ##   +/0 (type="-", u1 < z <= u2)
    cross1 <- vcross.lev[x$ev[,1]] ## left neighbour
    cross2 <- vcross.lev[x$ev[,2]] ## right neighbour
    vv.edges <- which((cross1 > 0) & (cross2 > 0) & (cross1 == cross2))
    for (edge in vv.edges) {
        lev <- cross1[edge]
        v1 <- x$tv[x$et[edge],x$eti[edge]]
        sign1 <- (z[v1] > levels[lev]+tol) - (z[v1] < levels[lev]-tol)
        neighb.t <- x$tt[x$et[edge],x$eti[edge]]
        if (is.na(neighb.t)) {
            v2 <- NA
            sign2 <- NA
        } else {
            v2 <- x$tv[neighb.t, x$tti[x$et[edge],x$eti[edge]] ]
            sign2 <- (z[v2] > levels[lev]+tol) - (z[v2] < levels[lev]-tol)
        }
        if (is.na(neighb.t)) {
            if (sign1 == 0) {
                idx <- rbind(idx, x$ev[edge,])
                grp <- c(grp, lev*2L + (type=="+")-(type=="-"))
            }
        } else if (((sign1 > 0) && (sign2 < 0)) ||
                   ((type=="+") && ((sign1 == 0) && (sign2 < 0))) ||
                   ((type=="-") && ((sign1 > 0) && (sign2 == 0))) ) {
            idx <- rbind(idx, x$ev[edge,])
            grp <- c(grp, lev*2L)
        }
    }

    ## For each boundary edge entirely between levels
    ##   store edge
    e.lower <- ecross.grp.lower + ((ecross.grp.lower-1) %% 2)
    e.upper <- ecross.grp.upper - ((ecross.grp.upper-1) %% 2)
    e.on.bnd <- which(is.na(x$tti[x$et+(x$eti-1)*x$Nt]))
    e.noncrossing <- e.on.bnd[ e.lower[e.on.bnd] == e.upper[e.on.bnd] ]
    idx <- rbind(idx, x$ev[e.noncrossing,,drop=FALSE])
    grp <- c(grp, as.integer(e.lower[e.noncrossing]))

    ## For each edge crossing at least one level,
    ##   calculate splitting vertices
    ##   if boundary edge, store new split edges
    e.crossing <- which(e.lower < e.upper)
    e.bndcrossing <- e.on.bnd[ e.lower[e.on.bnd] < e.upper[e.on.bnd] ]
    n$newv <- (sum(e.upper[e.crossing]-e.lower[e.crossing])/2 +
               sum(e.upper[e.bndcrossing]-e.lower[e.bndcrossing])/2)/2
    loc.new <- matrix(NA, n$newv, ncol(loc))
    loc.last <- 0L
    e.newv <- sparseMatrix(i=integer(0), j=integer(0), x=double(0),
                           dims=c(x$Ne, length(levels)))
    for (edge in e.crossing) {
        is.boundary.edge <- is.na(x$tt[x$et[edge], x$eti[edge]])
        if (is.boundary.edge) {
            edge.reverse <- NA
        } else {
        ## Interior edges appear twice; handle each only once.
            if (x$ev[edge,1] > x$ev[edge,2]) {
                next
            }
            edge.reverse <- x$te[x$tt[x$et[edge], x$eti[edge]],
                                 x$tti[x$et[edge], x$eti[edge]]]
        }
        ## lev = (1-beta_lev) * z1 + beta_lev * z2
        ##     = z1 + beta_lev * (z2-z1)
        ## beta_lev = (lev-z1)/(z2-z1)
        e.levels <- ((e.lower[edge]+1)/2):((e.upper[edge]-1)/2)
        beta <- ((levels[e.levels] - z[x$ev[edge,1]]) /
                 (z[x$ev[edge,2]] - z[x$ev[edge,1]]))
        loc.new.idx <- loc.last+seq_along(e.levels)
        loc.new[loc.new.idx,] <-
            (as.matrix(1-beta) %*% loc[x$ev[edge,1],,drop=FALSE]+
             as.matrix(beta) %*% loc[x$ev[edge,2],,drop=FALSE])
        e.newv[edge, e.levels] <- Nv + loc.new.idx
        if (!is.boundary.edge) {
            e.newv[edge.reverse, e.levels] <- Nv + loc.new.idx
        } else { ## edge is a boundary edge, handle now
            ev <- x$ev[edge,]
            the.levels <- which(e.newv[edge,] > 0)
            the.loc.idx <- e.newv[edge, the.levels ]
            the.levels <- c(min(the.levels)-1L, the.levels)*2L+1L
            if (z[ev[1]] > z[ev[2]]) {
                the.loc.idx <- the.loc.idx[length(the.loc.idx):1]
                the.levels <- the.levels[length(the.levels):1]
            }

            idx <- rbind(idx, cbind(c(ev[1], the.loc.idx),
                                    c(the.loc.idx, ev[2])))
            grp <- c(grp, the.levels)
        }
        loc.last <- loc.last + length(e.levels)
    }
    loc <- rbind(loc, loc.new)
    Nv <- nrow(loc)

    ## For each triangle, find non-level edge crossings, and
    ##   store new vertex-edge crossing edges
    ##   store new edge-edge crossing edges
    tris <- unique(x$et[e.crossing])
    for (tri in tris) {
        ## connect vertex-edge
        v.lev <- vcross.lev[x$tv[tri,]]
        for (vi in which(v.lev > 0)) {
            opposite.edge <- x$te[tri,vi]
            opposite.v <- e.newv[opposite.edge, v.lev[vi]]
            if (opposite.v > 0) {
                ## v2 on the right, v3 on the left
                v123 <- x$tv[tri, ((vi+(0:2)-1) %% 3) + 1]
                if (z[v123[3]] > z[v123[1]]) {
                    idx <- rbind(idx, cbind(v123[1], opposite.v))
                } else {
                    idx <- rbind(idx, cbind(opposite.v, v123[1]))
                }
                grp <- c(grp, v.lev[vi]*2L)
            }
        }
        ## connect edge-edge
        for (ei in 1:3) {
            edge <- x$te[tri, ei]
            next.edge <- x$te[tri, (ei %% 3L) + 1L]
            e.lev <- which(e.newv[edge, ] > 0)
            e.lev <- e.lev[ e.newv[next.edge, e.lev] > 0 ]
            if (length(e.lev) > 0) {
                ## v1 on the left, v2 on the right
                v12 <- x$ev[edge,]
                if (z[ v12[1] ] > z[ v12[2] ]) {
                    idx <- rbind(idx, cbind(e.newv[edge, e.lev],
                                            e.newv[next.edge, e.lev]))
                } else {
                    idx <- rbind(idx, cbind(e.newv[next.edge, e.lev],
                                            e.newv[edge, e.lev]))
                }
                grp <- c(grp, e.lev*2L)
            }
        }
    }

    ## Filter out unused nodes
    reo <- sort(unique(as.vector(idx)))
    loc <- loc[reo,,drop=FALSE]
    ireo <- integer(Nv)
    ireo[reo] <- seq_along(reo)
    idx <- matrix(ireo[idx], nrow(idx), 2)

    list(loc=loc, idx=idx, grp=grp)
}







tricontourmap <- function(x, z, nlevels = 10,
                          levels = pretty(range(z, na.rm = TRUE), nlevels),
                          ...)
{
    message("tricontourmap")
    UseMethod("tricontourmap")
}

tricontourmap.inla.mesh <-
    function(x, z, nlevels = 10,
             levels = pretty(range(z, na.rm = TRUE), nlevels),
             ...)
{
    message("tricontourmap.inla.mesh")
    tricontourmap.list(x$graph, z=z,
                       nlevels=nlevels, levels=levels,
                       loc=x$loc, ...)
}

tricontourmap.matrix <-
    function(x, z, nlevels = 10,
             levels = pretty(range(z, na.rm = TRUE), nlevels),
             loc, ...)
{
    message("tricontour.matrix")
    tricontourmap.list(list(tv=x), z=z,
                       nlevels=nlevels, levels=levels,
                       loc=loc, ...)
}



## Returns val=list(loc, idx, grp), where
##   grp = 1,...,nlevels*2+1, level groups are even, 2,4,...
## Suitable for
##   inla.mesh.segment(val$loc, val$idx[val$grp==k], val$idx[val$grp==k])
##     (supports R2 and S2)
## and, for odd k=1,3,...,nlevels*2-1,nlevels*2+1,
##   seg <- outline.to.inla.mesh.segment(val, grp.ccw=c(k-1,k), grp.cw=c(k+1))
##   sp <- outline.to.sp(val, grp.ccw=c(k-1,k), grp.cw=c(k+1), ccw=FALSE)


tricontourmap.list <-
    function(x, z, nlevels = 10,
             levels = pretty(range(z, na.rm = TRUE), nlevels),
             loc, type=c("+", "-"), tol=1e-7,
             output=c("sp", "inla.mesh.segment"), ...)
{
    message("tricontour.list")
    type <- match.arg(type)
    output <- match.arg(output)
    nlevels <- length(levels)

    tric <- tricontour(x=x, z=z, nlevels=nlevels, levels=levels,
                       loc=loc, type=type, tol=tol, ...)

    out <- list(map=list(), contour=list())
    for (k in seq_len(nlevels+1L)*2L-1L) {
        if (output == "sp") {
            out$map <- c(out$map,
                         list(try(outline.to.sp(tric,
                                                grp.ccw=c(k-1,k),
                                                grp.cw=c(k+1),
                                                ccw=FALSE,
                                                closed=TRUE))))
        } else {
            out$map <- c(out$map,
                         list(outline.to.inla.mesh.segment(tric,
                                                           grp.ccw=c(k-1,k),
                                                           grp.cw=c(k+1))))
        }
    }
    for (k in seq_len(nlevels)) {
        if (output == "sp") {
            out$contour <- c(out$contour,
                             list(try(outline.to.sp(tric,
                                                    grp.ccw=k*2L,
                                                    ccw=FALSE,
                                                    closed=FALSE))))
        } else {
            out$contour <-
                c(out$contour,
                  list(try(outline.to.inla.mesh.segment(tric,
                                                        grp.ccw=k*2L))))
        }
    }

    out
}




## excursions --> E,F
## contourmap --> E,P0123,F
## simconf
## continterp(excurobj, grid or mesh, outputgrid(opt), alpha, method)
## gaussint
