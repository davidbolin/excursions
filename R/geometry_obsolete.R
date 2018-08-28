## geometry_obsolete.R
##
##   Copyright (C) 2014-2016 Finn Lindgren, David Bolin
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

## In-filled points at transitions should have G[i] == -1
## to get a conservative approximation.
## To get only an "over/under set", use a constant non-negative integer G
##   and let calc.complement=FALSE
probabilitymap.old <-
  function(mesh, F, level, G,
           calc.complement=TRUE,
           tol=1e-7,
           output=c("sp", "inla.mesh.segment"),
           method, ...)
{
  output <- match.arg(output)

  if (output == "sp") {
    if (!requireNamespace("sp", quietly=TRUE)) {
      stop("The 'sp' package is needed.")
    }
  } else {
    if (!requireNamespace("INLA", quietly=TRUE)) {
      stop("The 'INLA' package is needed.")
    }
  }

  spout <- list()
  inlaout <- list()

  ## Find individual avoidance/between-level/under/over sets.
  for (k in sort(unique(G[G >= 0]))) {
    active.triangles <-
      which(rowSums(matrix(G[mesh$graph$tv] == k, nrow(mesh$graph$tv), 3)) == 3)
    active.nodes.idx <- unique(as.vector(mesh$graph$tv[active.triangles,]))

    if (length(active.nodes.idx) >= 3) { ## Non-empty mesh subset
      active.nodes <- logical(nrow(mesh$loc))
      active.nodes[active.nodes.idx] <- TRUE
      submesh <- submesh.mesh(active.nodes, mesh)

      subF <- rep(NA, nrow(submesh$loc))
      subF[submesh$idx$loc[active.nodes]] <- F[active.nodes]

#      if (FALSE) { ## Debugging plots
#        op <- par(mfrow=c(2,1))
#        on.exit(par(op))

#        class(mesh) <- "inla.mesh"
#        mesh$n <- nrow(mesh$loc)
#        mesh$manifold <- "R2"
#        proj <- inla.mesh.projector(mesh)
#        image(proj$x, proj$y, inla.mesh.project(proj, field=exp(F)),
#              zlim=range(exp(F)))
#        if (length(spout) > 0) {
#          plot(sp::SpatialPolygons(spout), add=TRUE, col="blue")
#        }
#        proj <- inla.mesh.projector(submesh)
#        image.plot(proj$x, proj$y, inla.mesh.project(proj, field=exp(subF)),
#                   xlim=range(mesh$loc[,1]), ylim=range(mesh$loc[,1]),
#                   zlim=range(exp(F)))
#        plot(submesh, add=TRUE)
#      }

      if (method == "step") {
        tric <- tricontour_step(x=submesh$graph, z=subF, levels=level,
                                loc=submesh$loc)
      } else {
        tric <- tricontour(x=submesh$graph, z=subF, levels=level,
                           loc=submesh$loc, type="+", tol=tol, ...)
      }
      ID <- as.character(k)

      if (output == "sp" || calc.complement) {
        spobj <- tryCatch(as.sp.outline(tric,
                                        grp.ccw=c(2,3),
                                        grp.cw=c(),
                                        ccw=FALSE,
                                        closed=TRUE,
                                        ID=ID),
                          error=function(e) NULL)
        if (!is.null(spobj)) {
          if (spobj@area == 0) {
            warning("Skipping zero area polygon in probabilitymap.")
          } else {
            spout[[ID]] <- spobj
          }
        }
      }
      if (output == "inla.mesh.segment") {
        inlaout[[ID]] <-
          as.inla.mesh.segment.outline(tric,
                                       grp.ccw=c(2,3),
                                       grp.cw=c(),
                                       grp=k)
      }
    }
}

  if (calc.complement) {
    ## Find contour set
    if (!requireNamespace("rgeos", quietly=TRUE)) {
      stop("Package 'rgeos' required for set complement calculations.")
    }

    ID <- "-1"
    outline <- INLA::inla.mesh.boundary(mesh)[[1]]
    sp.domain <- as.sp.outline(outline,
                               grp.ccw=unique(outline$grp),
                               grp.cw=integer(0),
                               ID=ID,
                               closed=TRUE)
    sp.domain <- sp::SpatialPolygons(list(sp.domain))

    if (length(spout) == 0) {
      ## Complement is the entire domain
      spout[[ID]] <- sp.domain@polygons[[1]]
    } else {
      spout.joined <- sp::SpatialPolygons(spout)
      spout.union <- rgeos::gUnaryUnion(spout.joined)
      spout[[ID]] <- rgeos::gDifference(sp.domain, spout.union)
      spout[[ID]] <- spout[[ID]]@polygons[[1]]
    }
    spout[[ID]]@ID <- ID

    if (output == "inla.mesh.segment") {
      inlaout[[ID]] <- INLA::inla.sp2segment(spout[[ID]])
    }
  }

  if (length(spout) > 0) {
    if (output == "sp") {
      out <- sp::SpatialPolygons(spout)
    } else {
      out <- do.call(INLA::inla.mesh.segment, inlaout)
    }
  } else {
    out <- NULL
  }

  out
}







continuous.old <- function(ex,
                       geometry,
                       alpha,
                       method=c("log", "logit", "linear", "step"),
                       output=c("sp", "inla"),
                       subdivisions=1,
                       calc.credible=TRUE)
{
  stopifnot(inherits(ex, "excurobj"))
  method <- match.arg(method)
  output <- match.arg(output)

  if (!(ex$meta$calculation %in% c("excursions",
                                   "contourmap"))) {
    stop(paste("Unsupported calculation '",
               ex$meta$calculation, "'.", sep=""))
  }

  if (missing(alpha)) {
    alpha <- ex$meta$alpha
  }
  if (alpha > ex$meta$F.limit) {
    warning(paste("Insufficient data: alpha = ", alpha,
                  " > F.limit = ", ex$meta$F.limit, sep=""))
  }

  info <- get.geometry(geometry)
  if (!(info$manifold %in% c("R2"))) {
    stop(paste("Unsupported manifold type '", info$manifold, "'.", sep=""))
  }
  if (length(ex$F) != prod(info$dims)) {
    stop(paste("The number of computed F-values (", length(ex$F), ") must match \n",
               "the number of elements of the continuous geometry definition (",
               prod(info$dims), ").", sep=""))
  }

  if (ex$meta$type == "=") {
    type <- "!="
    F.ex <- 1-ex$F
  } else {
    type <- ex$meta$type
    F.ex <- ex$F
  }
  F.ex[is.na(F.ex)] <- 0

  if (is.null(ex$meta$ind)) {
    active.nodes <- rep(TRUE, length(ex$F))
  } else {
    active.nodes <- logical(length(ex$F))
    active.nodes[ex$meta$ind] <- TRUE
  }
  if (info$geometry == "mesh") {
    mesh <- submesh.mesh(active.nodes, geometry)
  } else if (info$geometry == "lattice") {
    mesh <- submesh.grid(active.nodes, geometry)
  }
  mesh$graph <-
    generate.trigraph.properties(mesh$graph, Nv=nrow(mesh$loc))

  active.nodes <- !is.na(mesh$idx$loc)
  F.ex[mesh$idx$loc[active.nodes]] <- F.ex[active.nodes]
  G.ex <- rep(-1, nrow(mesh$loc))
  G.ex[mesh$idx$loc[active.nodes]] <- ex$G[active.nodes]

  ## Construct interpolation mesh
  F.geometry <- mesh
  F.geometry.A <- list()
  for (subdivision in seq_len(subdivisions)) {
    F.geometry <- subdivide.mesh(F.geometry)
    F.geometry.A <- c(F.geometry.A, list(F.geometry$A))
  }

  if (method == "log") {
    F.zero <- -1e20
    F.ex <- log(F.ex)
    level <- log(1-alpha)
    F.ex[is.infinite(F.ex) & F.ex < 0] <- F.zero
  } else if (method == "logit") {
    F.zero <- -1e20
    F.one <- +1e20
    F.ex <- log(F.ex)-log(1-F.ex)
    level <- log(1-alpha)-log(alpha)
    F.ex[is.infinite(F.ex) & F.ex < 0] <- F.zero
    F.ex[is.infinite(F.ex) & F.ex > 0] <- F.one
  } else if (method == "linear") {
    F.zero <- 0
    level <- 1-alpha
  } else {
    ## 'step'
    F.zero <- 0
    level <- 1-alpha
  }

  ## For ordinary excursions, the input set/group information is not used.
  if (type == ">") {
    G.ex <- rep(1, length(G.ex))
  } else if (type == "<") {
    G.ex <- rep(0, length(G.ex))
  }
  ## Copy 'G' and interpolate 'F' within coherent single level regions.
  G.interp <- G.ex
  F.interp <- F.ex
  for (subdivision in seq_len(subdivisions)) {
    G.input <- G.interp
    F.input <- F.interp
    G.interp <- rep(-1, nrow(F.geometry.A[[subdivision]]))
    F.interp <- rep(F.zero, nrow(F.geometry.A[[subdivision]]))
    for (k in unique(G.ex[G.ex >= 0])) {
      ok.in <- (G.input == k)
      ok.out <- (rowSums(F.geometry.A[[subdivision]][,ok.in,drop=FALSE]) >
                 1 - 1e-12)
      G.interp[ok.out] <- k
    }
    ok.in <- (G.input >= 0)
    ok.out <- (G.interp >=0)
    if (method =="step") {
      for (vtx in which(ok.out)) {
        F.interp[vtx] <- min(F.input[F.geometry.A[[subdivision]][vtx, ] > 0])
      }
    } else {
      F.interp[ok.out] <-
        as.vector(F.geometry.A[[subdivision]][ok.out,ok.in,drop=FALSE] %*%
                  F.input[ok.in])
    }
  }

  if (method == "log") {
    F.interp.nontransformed <- exp(F.interp)
  } else if (method == "logit") {
    F.interp.nontransformed <- 1/(1 + exp(-F.interp))
  } else if (method == "linear") {
    F.interp.nontransformed <- F.interp
  } else {
    ## 'step'
    F.interp.nontransformed <- F.interp
  }
  F.interp.nontransformed[G.interp == -1] <- 0
  if (ex$meta$type == "=") {
    F.interp.nontransformed <- 1-F.interp.nontransformed
  }

  M <- probabilitymap.old(F.geometry,
                          F=F.interp,
                          level=level,
                          G=G.interp,
                          calc.complement=calc.credible,
                          method=method,
                          output=output)

  if (requireNamespace("INLA", quietly=TRUE)) {
    F.geometry <- INLA::inla.mesh.create(loc=F.geometry$loc,
                                         tv=F.geometry$graph$tv)
    ## Handle possible node reordering in inla.mesh.create()
    F.interp.nontransformed[F.geometry$idx$loc] <- F.interp.nontransformed
    G.interp[F.geometry$idx$loc] <- G.interp
  }

  out <- list(F=F.interp.nontransformed, G=G.interp, M=M, F.geometry=F.geometry)

  if (!is.null(ex$P0)) {
    if (!requireNamespace("INLA", quietly=TRUE)) {
      warning("The 'INLA' package is required for P0 calculations.")
    } else {
      fem <- INLA::inla.mesh.fem(F.geometry, order=1)
      out$P0 <-
        sum(diag(fem$c0) * F.interp.nontransformed) /
        sum(diag(fem$c0))
    }
  }

  out
}
