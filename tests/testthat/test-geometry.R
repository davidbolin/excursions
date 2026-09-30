test_that("Contours on triangulations are unchanged", {
  tm <- testdata.mesh()
  tol <- 1e-10
  expect_equal(
    unlist(fingerprint.tricontour(
      tricontour(tm$mesh, z = tm$z, levels = c(-1, -0.3, 0.4, 1.1))
    )),
    REF$tc.plus,
    tolerance = tol
  )
  expect_equal(
    unlist(fingerprint.tricontour(
      tricontour(tm$mesh, z = tm$z, levels = c(-1, 0, 1), type = "-")
    )),
    REF$tc.minus,
    tolerance = tol
  )
  ## Many vertices exactly on the levels
  zl <- round(tm$z * 2) / 2
  expect_equal(
    unlist(fingerprint.tricontour(
      tricontour(tm$mesh, z = zl, levels = c(-0.5, 0, 0.5))
    )),
    REF$tc.onlevel,
    tolerance = tol
  )
})

test_that("Contour edges are consistent", {
  tm <- testdata.mesh(15)
  levels <- c(-1, 0, 1)
  tc <- tricontour(tm$mesh, z = tm$z, levels = levels)
  expect_equal(ncol(tc$idx), 2)
  expect_equal(nrow(tc$idx), length(tc$grp))
  expect_true(all(tc$idx >= 1 & tc$idx <= nrow(tc$loc)))
  expect_true(all(tc$grp >= 1 & tc$grp <= 2 * length(levels) + 1))
  ## Edges on a level have both end points on that level. The field is
  ## piecewise linear, so interpolate it with the mesh basis.
  A <- fmesher::fm_basis(tm$mesh, loc = tc$loc)
  z.loc <- as.vector(A %*% tm$z)
  on.level <- tc$grp %% 2 == 0
  lev <- levels[tc$grp[on.level] / 2]
  expect_equal(z.loc[tc$idx[on.level, 1]], lev, tolerance = 1e-8)
  expect_equal(z.loc[tc$idx[on.level, 2]], lev, tolerance = 1e-8)
})

## Check that the sequences from connect.segments follow the segments
check.segments <- function(res, segments, grp, grp.ccw, grp.cw, ccw = TRUE) {
  oriented <- segments
  cw <- grp %in% grp.cw
  oriented[cw, ] <- segments[cw, 2:1]
  used <- unlist(res$seg)
  expect_setequal(used, which(grp %in% c(grp.ccw, grp.cw)))
  expect_false(anyDuplicated(used) > 0)
  for (k in seq_along(res$sequences)) {
    s <- res$sequences[[k]]
    seg <- res$seg[[k]]
    if (!ccw) {
      s <- rev(s)
      seg <- rev(seg)
    }
    expect_equal(length(s), length(seg) + 1)
    expect_equal(oriented[seg, 1], s[-length(s)])
    expect_equal(oriented[seg, 2], s[-1])
    expect_equal(res$grp[[k]], if (ccw) grp[seg] else rev(grp[seg]))
  }
}

test_that("Segments are connected into sequences", {
  ## A closed loop, an open chain where the first segment is in the middle,
  ## and a chain with two segments
  seg <- rbind(
    c(1, 2), c(2, 3), c(3, 4), c(4, 1),
    c(5, 6), c(7, 5), c(6, 8),
    c(10, 11), c(9, 10)
  )
  grp <- c(1, 1, 2, 2, 1, 1, 1, 3, 3)

  res <- excursions:::connect.segments(seg, segment.grp = grp)
  expect_equal(unlist(fingerprint.segments(res)), REF$seg1)
  check.segments(res, seg, grp, grp.ccw = unique(grp), grp.cw = integer(0))
  expect_equal(res$sequences[[1]], c(1, 2, 3, 4, 1))
  expect_equal(res$sequences[[2]], c(7, 5, 6, 8))
  expect_equal(res$sequences[[3]], c(9, 10, 11))

  res <- excursions:::connect.segments(seg, grp,
    grp.ccw = c(1, 3), grp.cw = 2, ccw = FALSE
  )
  expect_equal(unlist(fingerprint.segments(res)), REF$seg2)
  check.segments(res, seg, grp, grp.ccw = c(1, 3), grp.cw = 2, ccw = FALSE)

  ## Contour segments on a mesh
  tm <- testdata.mesh()
  tc <- tricontour(tm$mesh, z = tm$z, levels = c(-1, 0, 1))
  res <- excursions:::connect.segments(tc$idx, tc$grp, grp.ccw = c(2, 3), grp.cw = 4)
  expect_equal(unlist(fingerprint.segments(res)), REF$seg3)
  check.segments(res, tc$idx, tc$grp, grp.ccw = c(2, 3), grp.cw = 4)

  ## No segments
  res <- excursions:::connect.segments(seg, grp, grp.ccw = 5)
  expect_length(res$sequences, 0)
})

test_that("Continuous excursion sets are unchanged", {
  skip_if_not_installed("sp")
  d <- testdata.spde(10)
  lat <- fmesher::fm_lattice_2d(x = d$x, y = d$x)
  ex <- excursions(
    alpha = 0.1, u = 0.5, mu = d$mu, Q = d$Q, type = ">",
    seed = d$seed, max.threads = 1, n.iter = 2000
  )
  for (method in c("step", "linear", "log")) {
    r <- continuous(ex, lat, method = method, output = "sp")
    crd <- unlist(lapply(r$M@polygons, function(p) {
      lapply(p@Polygons, function(q) q@coords)
    }))
    expect_equal(
      c(
        F = ref.summary(r$F), P0 = r$P0, M = length(crd), Msum = sum(crd),
        Mw = sum(crd * seq_along(crd))
      ),
      REF[[paste0("cont.", method)]],
      tolerance = 1e-10
    )
  }
})
