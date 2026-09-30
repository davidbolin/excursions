test_that("Variances", {
  data <- integration.testdata1()
  vars <- excursions.variances(data$L, max.threads = 1)
  v <- diag(solve(data$Q))
  expect_equal(vars, v, tolerance = 1e-7)
})


test_that("Variances Q and L", {
  data <- integration.testdata1()
  v1 <- excursions.variances(L = data$L, max.threads = 1)
  v2 <- excursions.variances(Q = data$Q, max.threads = 1)
  expect_equal(v1, v2, tolerance = 1e-7)
})


test_that("Variances for the different forms of the Cholesky factor", {
  data <- testdata.spde(15)
  Q <- data$Q
  v <- diag(solve(as.matrix(Q)))

  ## Precision matrix, factorised with a fill reducing permutation
  expect_equal(excursions.variances(Q = Q), v, tolerance = 1e-10)

  ## Lower and upper triangular sparse factors without permutation
  L <- as(Matrix::Cholesky(Q, perm = FALSE, LDL = FALSE), "CsparseMatrix")
  expect_equal(L@uplo, "L")
  expect_equal(excursions.variances(L = L), v, tolerance = 1e-10)
  expect_equal(excursions.variances(L = t(L)), v, tolerance = 1e-10)
  expect_equal(excursions.variances(L = chol(Q)), v, tolerance = 1e-10)

  ## Supernodal factor, which may store structural zeros
  Ls <- as(Matrix::Cholesky(Q, perm = FALSE, LDL = FALSE, super = TRUE), "CsparseMatrix")
  expect_equal(excursions.variances(L = Ls), v, tolerance = 1e-10)

  ## Dense factor
  expect_equal(excursions.variances(L = chol(as.matrix(Q))), v, tolerance = 1e-10)

  ## A permuted factor gives the variances in the order of the factor
  ch <- Matrix::Cholesky(Q, perm = TRUE, LDL = FALSE)
  Lp <- as(ch, "CsparseMatrix")
  perm <- ch@perm + 1L
  expect_false(identical(perm, seq_len(data$n)))
  expect_equal(excursions.variances(L = Lp), v[perm], tolerance = 1e-10)

  ## The number of threads does not change the result
  expect_identical(
    excursions.variances(L = L, max.threads = 1),
    excursions.variances(L = L, max.threads = 4)
  )
})

test_that("Variances for special cases", {
  ## Unit diagonal triangular factor
  L <- methods::new("dtCMatrix",
    i = c(1L, 2L, 2L), p = c(0L, 2L, 3L, 3L), x = c(0.5, -0.2, 0.3),
    Dim = c(3L, 3L), uplo = "L", diag = "U"
  )
  v <- diag(solve(as.matrix(L %*% t(L))))
  expect_equal(excursions.variances(L = L), v, tolerance = 1e-12)

  ## One node
  expect_equal(excursions.variances(Q = Matrix::Matrix(4, 1, 1, sparse = TRUE)), 0.25)
  expect_equal(excursions.variances(L = Matrix::Matrix(2, 1, 1, sparse = TRUE)), 0.25)

  ## Diagonal precision
  Q <- Matrix::Diagonal(5, 1:5)
  expect_equal(excursions.variances(Q = Q), 1 / (1:5))
  expect_equal(excursions.variances(L = sqrt(Q)), 1 / (1:5))
})

test_that("Variances give an error for a factor without a diagonal element", {
  L <- Matrix::sparseMatrix(
    i = c(1, 2, 3), j = c(1, 1, 3), x = c(1, 0.5, 1),
    dims = c(3, 3), triangular = TRUE
  )
  expect_error(excursions.variances(L = L), "diagonal")
})
