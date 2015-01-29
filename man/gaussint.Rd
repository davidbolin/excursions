\name{gaussint}
\alias{gaussint}
\title{
  Sequential estimation of Gaussian integrals
}
\description{
  \code{gaussint} is used for calculating Gaussian integrals
 \deqn{latex}{
  \int_a^b \frac{|Q|^{1/2}}{(2\pi)^{n/2}} \exp(-\frac1{2}(x-\mu)^{T}Q(x-\mu)) dx
 }
 A limit value \eqn{lim} can be used to stop the integration if the
 sequential estimate goes below the limit, which can result in substantial
 computational savings in cases when one only is interested in testing if the
 integral is above the limit value. The integral is calculated sequentially, and
 estimates for all subintegrals are also returned.

}
\usage{
gaussint(mu,
         Q.chol,
         Q,
         a,
         b,
         lim = 0,
         n.iter = 10000,
         ind,
         use.reordering = c("natural","sparsity","limits"),
         max.size,
         max.threads=0,
         seed,
         LDL=FALSE)
}
\arguments{
  \item{mu }{Expectation vector for the Gaussian distribution.}
  \item{Q.chol }{The Cholesky factor of the precision matrix
	  (optional)}
  \item{Q }{Precision matrix for the Gaussian distribution. If Q is supplied but not Q.chol, the cholesky factor is computed before integrating.}
  \item{a }{Lower limit in integral.}
  \item{b }{Upper limit in integral.}
  \item{lim }{If this argument is used, the integration
    is stopped and 0 is returned if the estimated value goes
    below \eqn{lim}.}
  \item{n.iter }{Number or iterations in the MC sampler that is
    used for approximating probabilities. The default value is
    10000.}
  \item{ind }{Indices of the nodes that should be analyzed (optional)}
  \item{use.reordering }{Determines what reordering to use:
    \itemize{
      \item{"natural" }{No reordering is performed.}
      \item{"sparsity" }{Reorder for sparsity in the cholesky factor (MMD reordering is used).}
      \item{"limits" }{Reorder by moving all nodes with a=-Inf and b=Inf first and then reordering for sparsity (CAMD reordering is used).}
    }}
  \item{max.size }{The largest number of sub-integrals to compute. Default is the total dimension of the distribution. }
  \item{max.threads }{Decides the number of threads the program can use. Set to 0 for using the maximum number of threads allowed by the system (default).}
  \item{seed }{The random seed to use (optional). }
  \item{LDL }{Use LDL factorisations? This is useful for analysis of problems with positive semidefinite precisions.}

}
\value{
  A list:
  \item{P }{Value of the integral.}
  \item{E }{Estimated error of the P estimate.}
  \item{Pv }{A vector with the estimates of all sub-integrals.}
  \item{Ev }{A vector with the estimated errors of the Pv estimates}
}
\author{
  David Bolin \email{davidbolin@gmail.com}
}

\references{
Bolin, D. and Lindgren, F. (2013) \emph{Excursion and contour uncertainty regions for latent Gaussian models}, Journal of the Royal Statistical Society: Series B (in press).
}

\examples{
## Create mean and a tridiagonal precision matrix
n = 11
mu.x = seq(-5, 5, length=n)
Q.x = Matrix(toeplitz(c(1, -0.1, rep(0, n-2))))

## Calculate the probability that the process is between mu-3 and mu+3
prob = gaussint(mu=mu.x, Q=Q.x, a= mu.x-3, b=mu.x+3, max.threads=2)
prob$P
}
