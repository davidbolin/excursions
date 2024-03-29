% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gaussint.R
\name{gaussint}
\alias{gaussint}
\title{Sequential estimation of Gaussian integrals}
\usage{
gaussint(
  mu,
  Q.chol,
  Q,
  a,
  b,
  lim = 0,
  n.iter = 10000,
  ind,
  use.reordering = c("natural", "sparsity", "limits"),
  max.size,
  max.threads = 0,
  seed
)
}
\arguments{
\item{mu}{Expectation vector for the Gaussian distribution.}

\item{Q.chol}{The Cholesky factor of the precision matrix (optional).}

\item{Q}{Precision matrix for the Gaussian distribution. If Q is supplied but not Q.chol,
the cholesky factor is computed before integrating.}

\item{a}{Lower limit in integral.}

\item{b}{Upper limit in integral.}

\item{lim}{If this argument is used, the integration is stopped and 0 is returned
if the estimated value goes below \eqn{lim}.}

\item{n.iter}{Number or iterations in the MC sampler that is used for approximating
probabilities. The default value is 10000.}

\item{ind}{Indices of the nodes that should be analyzed (optional).}

\item{use.reordering}{Determines what reordering to use:
\describe{
  \item{"natural" }{No reordering is performed.}
  \item{"sparsity" }{Reorder for sparsity in the cholesky factor (MMD reordering
  is used).}
  \item{"limits" }{Reorder by moving all nodes with a=-Inf and b=Inf first and
  then reordering for sparsity (CAMD reordering is used).}
  }}

\item{max.size}{The largest number of sub-integrals to compute. Default is the total
dimension of the distribution.}

\item{max.threads}{Decides the number of threads the program can use. Set to 0 for
using the maximum number of threads allowed by the system (default).}

\item{seed}{The random seed to use (optional).}
}
\value{
A list with elements
\item{P }{Value of the integral.}
\item{E }{Estimated error of the P estimate.}
\item{Pv }{A vector with the estimates of all sub-integrals.}
\item{Ev }{A vector with the estimated errors of the Pv estimates.}
}
\description{
\code{gaussint} is used for calculating \eqn{n}-dimensional Gaussian integrals
\deqn{\int_a^b \frac{|Q|^{1/2}}{(2\pi)^{n/2}}
\exp(-\frac1{2}(x-\mu)^{T}Q(x-\mu)) dx}{|Q|^(1/2)*(2\pi)^(-n/2) \int_a^b exp(-0.5*(x-\mu)^T Q (x-\mu)) dx}
A limit value \eqn{lim} can be used to stop the integration if the sequential
estimate goes below the limit, which can result in substantial computational
savings in cases when one only is interested in testing if the integral is above
the limit value. The integral is calculated sequentially, and estimates for
all subintegrals are also returned.
}
\details{
The function uses sequential importance sampling to estimate the
Gaussian integral, and returns all computed sub-integrals. This means that if, for
example, the function is used to compute \eqn{P(x>0)} for an n-dimensional Gaussian
variable \eqn{x}, then all integrals \eqn{P(x_1>0,\ldots,x_i>0)} for \eqn{i=1,\ldots,n} are
computed.

If one is only interested in whether \eqn{P(x>0)>\alpha} or not, then one can
stop the integration as soon as \eqn{P(x_1>0,\ldots,x_i>0)<\alpha}. This can save a lot of
computation time if \eqn{P(x_1>0,\ldots,x_i>0)< \alpha} for \eqn{i} much smaller than
\eqn{n}. This limit value is specified by the \code{lim} argument.

Which reordering to use depends on what the purpose of the calculation is and what
the integration limits are. However, in general the \code{limits} reordering is typically
most appropriate since this combines sparisty (which improves accuracy and reduces
computational cost) with automatic handling of dimensions with limits \code{a=-Inf} and
\code{b=Inf}, which do not affect the probability but affect the computation time
if they are not handled separately.
}
\examples{
## Create mean and a tridiagonal precision matrix
n <- 11
mu.x <- seq(-5, 5, length = n)
Q.x <- Matrix(toeplitz(c(1, -0.1, rep(0, n - 2))))
## Calculate the probability that the variable is between mu-3 and mu+3
prob <- gaussint(mu = mu.x, Q = Q.x, a = mu.x - 3, b = mu.x + 3, max.threads = 2)
prob$P
}
\references{
Bolin, D. and Lindgren, F. (2015) \emph{Excursion and contour uncertainty regions for latent Gaussian models}, JRSS-series B, vol 77, no 1, pp 85-106.

Bolin, D. and Lindgren, F. (2018), \emph{Calculating Probabilistic Excursion Sets and Related Quantities Using excursions}, Journal of Statistical Software, vol 86, no 1, pp 1-20.
}
\author{
David Bolin \email{davidbolin@gmail.com}
}
