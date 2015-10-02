\name{simconf.mc}
\alias{simconf.mc}
\title{
  Simultaneous confidence regions using Monte Carlo samples
}
\description{
  \code{simconf.mc} is used for calculating simultaneous
  confidence regions based on Monte Carlo samples. The function returns
  upper and lower bounds \eqn{a} and \eqn{b} such that \eqn{P(a<x<b) =
    1-alpha}.
}
\usage{
simconf.mc(samples,
           alpha,
           ind,
           verbose=FALSE)
}
\arguments{
  \item{samples}{Matrix with model Monte Carlo samples. Each column contains a sample of the model.}
  \item{alpha}{Error probability for the region.}
  \item{ind}{Indices of the nodes that should be analyzed
    (optional)}
  \item{verbose}{Set to TRUE for verbose mode (optional)}
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
}

\examples{
## Create mean and a tridiagonal precision matrix
n = 11
mu.x = seq(-5, 5, length=n)
Q.x = Matrix(toeplitz(c(1, -0.1, rep(0, n-2))))

## Sample the model 100 times (increase for better estimate)
X = mu.x + solve(chol(Q.x),matrix(rnorm(n=n*100),nrow=n,ncol=100))

## calculate the confidence region
conf = simconf.mc(X,0.2)
## Plot the region
plot(mu.x, type="l", ylim=c(-10, 10),
     main='Mean (black) and confidence region (red)')
lines(conf$a, col=2)
lines(conf$b, col=2)
}
