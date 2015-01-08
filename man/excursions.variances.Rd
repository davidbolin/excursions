\name{excursions.variances}
\alias{excursions.variances}
\title{
  Calculate variances from a sparse precision matrix
}
\description{
  \code{excursions.variances} calculates the diagonal of the inverse of a sparse symmetric positive definite matrix \code{Q}. The method for calculating the diagonal requires the Cholesky factor of \code{Q}, which should be supplied if available.
}
\usage{
excursions.variances(L,Q)
}
\arguments{
  \item{L}{Cholesky factor of precision matrix.}
  \item{Q}{Precision matrix.}
}
\value{
  A vector with the diagonal elements.}
\author{
  David Bolin \email{davidbolin@gmail.com}
}

\examples{
## Create a tridiagonal precision matrix
n = 21
Q = Matrix(toeplitz(c(1, -0.1, rep(0, n-2))))

v2 = excursions.variances(Q=Q)

## var2 should be the same as:
v1 = diag(solve(Q))

}
