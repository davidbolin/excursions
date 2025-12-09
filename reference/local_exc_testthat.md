# Unit test helpers

Local helper functions for package unit tests

## Usage

``` r
local_exc_safe_inla(multicore = FALSE, quietly = TRUE, envir = parent.frame())
```

## Arguments

- multicore:

  logical; if `TRUE`, multiple cores are allowed, and the INLA
  `num.threads` option is not checked or altered. Default: `FALSE`,
  multicore not allowed (used for examples and unit tests).

- quietly:

  logical; if `TRUE`, prints diagnostic messages. A message is always
  printed if the INLA `num.threads` option is altered, regardless of the
  `quietly` argument. Default: TRUE.

- envir:

  environment for exit handlers

## Functions

- `local_exc_safe_inla()`: Tests should set num.threads = "1:1" to
  ensure within-system repeatability by calling `local_exc_safe_inla()`;
  see also
  [`exc_safe_inla()`](https://davidbolin.github.io/excursions/reference/exc_safe_inla.md)

## Examples

``` r
if (FALSE) { # \dontrun{
local_exc_safe_inla(multicore = FALSE)
} # }
```
