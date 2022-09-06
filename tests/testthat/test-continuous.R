test_that("Continous on contourmap, R2 mesh", {
  skip_on_cran()
  local_exc_safe_inla()

    data <- integration.testdata1()
    res1 = contourmap(data$mu, data$Q, n.levels=2,
                      seed=data$seed, alpha=0.1,max.threads=1,
                      compute = list(F = TRUE,
                                     measures = c("P2","P1","P0")))
    
    bnd <- INLA::inla.mesh.segment(cbind(c(0, 1, 2, 2, 2, 1, 0, 0),
                                         c(0, 0, 0, 1, 2, 2, 2, 1)),
                                   is.bnd = TRUE)
    loc <- cbind(c(0.5, 1.5, 1),
                 c(0.5, 0.5, 1.5))
    mesh <- INLA::inla.mesh.create(loc = loc, boundary = bnd)
    expect_equal(data$n, mesh$n)
    
    GQ_direct <- gaussquad(mesh, method = "direct")
    GQ_make_A <- gaussquad(mesh, method = "make.A")
    expect_equal(GQ_direct$A, GQ_make_A$A, tolerance = 1e-15)
    
    res2 <- continuous(res1, mesh, method = "linear")
    
    expect_error(continuous(res1, mesh, method = "linear", output = "inla"),
                 "Output format 'inla' not supported for 'calc.credible = TRUE'.")
    res3 <- continuous(res1, mesh, method = "linear", output = "inla",
                       calc.credible = FALSE)
    res4 <- continuous(res1, mesh, method = "log", output = "inla",
                       calc.credible = FALSE)
    res5 <- continuous(res1, mesh, method = "step", output = "inla",
                       calc.credible = FALSE)
    
    expect_s4_class(res2$M, "SpatialPolygons")
    expect_s3_class(res3$M, "inla.mesh.segment")
    expect_s3_class(res4$M, "inla.mesh.segment")
    expect_s3_class(res5$M, "inla.mesh.segment")
})

test_that("Continous on contourmap, M mesh", {
  skip_on_cran()
  local_exc_safe_inla()

    data <- integration.testdata1()
    res1 = contourmap(data$mu, data$Q, n.levels=2,
                      seed=data$seed, alpha=0.1,max.threads=1)
    
    bnd <- INLA::inla.mesh.segment(cbind(c(0, 1, 2, 2, 2, 1, 0, 0),
                                         c(0, 0, 0, 1, 2, 2, 2, 1)),
                                   is.bnd = TRUE)
    loc <- cbind(c(0.5, 1.5, 1),
                 c(0.5, 0.5, 1.5))
    mesh <- INLA::inla.mesh.create(loc = loc, boundary = bnd)
    expect_equal(data$n, mesh$n)
    
    # Alter z-coordinates and mark as general manifold
    mesh$loc[,3] = seq_len(mesh$n)
    mesh$manifold = "M"
    
    res2 <- continuous(res1, mesh, method = "linear",
                       output = "inla",
                       calc.credible = FALSE)
    
    expect_s3_class(res2$M, "inla.mesh.segment")
    
    res3 <- continuous(res1, mesh, method = "log",
                       output = "inla",
                       calc.credible = FALSE)
    
    expect_s3_class(res3$M, "inla.mesh.segment")
    
    res4 <- continuous(res1, mesh, method = "step",
                       output = "inla",
                       calc.credible = FALSE)
    
    expect_s3_class(res4$M, "inla.mesh.segment")

})
