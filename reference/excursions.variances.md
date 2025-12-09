# Calculate variances from a sparse precision matrix

`excursions.variances` calculates the diagonal of the inverse of a
sparse symmetric positive definite matrix `Q`.

## Usage

``` r
excursions.variances(L, Q, max.threads = 0)
```

## Arguments

- L:

  Cholesky factor of precision matrix.

- Q:

  Precision matrix.

- max.threads:

  Decides the number of threads the program can use. Set to 0 for using
  the maximum number of threads allowed by the system (default).

## Value

A vector with the variances.

## Details

The method for calculating the diagonal requires the Cholesky factor,
`L`, of `Q`, which should be supplied if available. If `Q` is provided,
the cholesky factor is calculated and the variances are then returned in
the same ordering as `Q`. If `L` is provided, the variances are returned
in the same ordering as `L`, even if `L@invpivot` exists.

## Author

David Bolin <davidbolin@gmail.com>

## Examples

``` r
## Create a tridiagonal precision matrix
n <- 21
Q <- Matrix(toeplitz(c(1, -0.1, rep(0, n - 2))))
v2 <- excursions.variances(Q = Q, max.threads = 2)
## var2 should be the same as:
v1 <- diag(solve(Q))
```
