\name{excursions.variances}
\alias{excursions.variances}
\title{
  Calculate variances from a sparse precision matrix
}
\description{
  \code{excursions.variances} calculates the diagonal of the inverse of a sparse symmetric positive definite matrix \code{Q}. The method for calculating the diagonal requires the Cholesky factor of \code{Q}, which should be supplied if available.
}
\usage{
excursions.variances(L,Q,max.threads=0,use.spam=FALSE)
}
\arguments{
  \item{L}{Cholesky factor of precision matrix.}
  \item{Q}{Precision matrix.}
  \item{max.threads}{Decides the number of threads the program can use. Set to 0 for using the maximum number of threads allowed by the system (default).}
  \item{use.spam}{If set to true TRUE and if a precision matrix is provided, SPAM is used for the Cholesky decomposition.}
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

v2 = excursions.variances(Q=Q,max.threads=2)

## var2 should be the same as:
v1 = diag(solve(Q))

}
