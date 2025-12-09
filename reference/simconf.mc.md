# Simultaneous confidence regions using Monte Carlo samples

`simconf.mc` is used for calculating simultaneous confidence regions
based on Monte Carlo samples. The function returns upper and lower
bounds \\a\\ and \\b\\ such that \\P(a\<x\<b) = 1-\alpha\\.

## Usage

``` r
simconf.mc(samples, alpha, ind, verbose = FALSE)
```

## Arguments

- samples:

  Matrix with model Monte Carlo samples. Each column contains a sample
  of the model.

- alpha:

  Error probability for the region.

- ind:

  Indices of the nodes that should be analyzed (optional).

- verbose:

  Set to TRUE for verbose mode (optional).

## Value

An object of class "excurobj" with elements

- a :

  The lower bound.

- b :

  The upper bound.

- a.marginal :

  The lower bound for pointwise confidence bands.

- b.marginal :

  The upper bound for pointwise confidence bands.

## Details

See
[`simconf()`](https://davidbolin.github.io/excursions/reference/simconf.md)
for details.

## See also

[`simconf()`](https://davidbolin.github.io/excursions/reference/simconf.md),
[`simconf.inla()`](https://davidbolin.github.io/excursions/reference/simconf.inla.md)

## Author

David Bolin <davidbolin@gmail.com>

## Examples

``` r
## Create mean and a tridiagonal precision matrix
n <- 11
mu.x <- seq(-5, 5, length = n)
Q.x <- Matrix(toeplitz(c(1, -0.1, rep(0, n - 2))))
## Sample the model 100 times (increase for better estimate)
X <- mu.x + solve(chol(Q.x), matrix(rnorm(n = n * 100), nrow = n, ncol = 100))
## calculate the confidence region
conf <- simconf.mc(X, 0.2)
#> 0.02017992
## Plot the region
plot(mu.x,
  type = "l", ylim = c(-10, 10),
  main = "Mean (black) and confidence region (red)"
)
lines(conf$a, col = 2)
lines(conf$b, col = 2)
```
