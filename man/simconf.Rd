\name{simconf}
\alias{simconf}
\title{
  Simultaneous confidence regions for Gaussian models
}
\description{
  \code{simconf} is used for calculating simultaneous
  confidence regions for Gaussian models \eqn{x}. The function returns
  upper and lower bounds \eqn{a} and \eqn{b} such that \eqn{P(a<x<b) =
    1-alpha}.
}
\usage{
simconf(alpha,
        mu,
        Q,
        n.iter=10000,
        Q.chol,
        vars,
        ind=NULL,
        verbose=0,
        max.threads=0,
        seed=NULL,
        LDL=TRUE)
}
\arguments{
  \item{alpha}{Error probability for the region.}
  \item{mu}{Expectation vector for the Gaussian distribution.}
  \item{Q}{Precision matrix for the Gaussian distribution.}
  \item{n.iter}{Number or iterations in the MC sampler that is
    used for approximating probabilities. The default value is
    10000.}
  \item{Q.chol}{The Cholesky factor of the precision matrix
    (optional)}
  \item{vars}{Precomputed marginal variances (optional)}
  \item{ind}{Indices of the nodes that should be analyzed
    (optional)}
  \item{verbose}{Set to TRUE for verbose mode (optional)}
  \item{max.threads}{Decides the number of threads the program can use. Set to 0 for using the maximum number of threads allowed by the system (default).}
  \item{seed}{Random seed (optional)}
  \item{LDL}{Use LDL factorisations? This is useful for analysis of problems with positive semidefinite precisions.}

}
\value{
  A list:
  \item{a }{The lower bound.}
  \item{b }{The upper bound.}
  \item{a.marginal }{The lower bound for pointwise confidence bands.}
  \item{b.marginal }{The upper bound for pointwise confidence bands.}

}
\author{
  David Bolin \email{davidbolin@gmail.com} and
  Finn Lindgren \email{finn.lindgren@gmail.com}
}

\references{
Bolin et al. Statistical prediction of global sea level from global temperature, Statistica Sinica, Vol 25, pp 351-367.
}

\examples{
## Create mean and a tridiagonal precision matrix
n = 11
mu.x = seq(-5, 5, length=n)
Q.x = Matrix(toeplitz(c(1, -0.1, rep(0, n-2))))

## calculate the confidence region
conf = simconf(0.05, mu.x, Q.x, max.threads=2)

## Plot the region
plot(mu.x, type="l", ylim=c(-10, 10),
     main='Mean (black) and confidence region (red)')
lines(conf$a, col=2)
lines(conf$b, col=2)
}
