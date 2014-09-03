\name{excursions.int}
\alias{excursions.int}
\title{
  Calculate Gaussian integrals
}
\description{
  \code{excursions.int} is used for calculating Gaussian integrals
 \deqn{latex}{
 	\int_a^b
	\frac{|Q|^{1/2}}{(2\pi)^{n/2}}
	\exp(-\frac1{2}(x-\mu)^{T}Q(x-\mu))
	dx
 } 
 A limit value \eqn{alpha} can be used to stop the integration if the
 value is below \eqn{1-alpha}.
}
\usage{
excursions.int(mu, Q, a, b, n.iter=10000, Q.chol, ind,
               verbose=0, alpha=1, max.threads=0)
}
\arguments{
  	\item{mu}{Expectation vector for the Gaussian distribution.}
  	\item{Q}{Precision matrix for the Gaussian distribution.}
  	\item{a}{Lower limit in integral.}
  	\item{b}{Upper limit in integral.}
  	\item{n.iter}{Number or iterations in the MC sampler that is
	  used for approximating probabilities. The default value is
	  10000.}
  	\item{Q.chol}{The Cholesky factor of the precision matrix
	  (optional)}
	\item{ind}{Indices of the nodes that should be analyzed
	  (optional)}
	\item{verbose}{Set to TRUE for verbose mode (optional)}
	\item{alpha}{If this argument is used, the integration
	  is stopped and 0 is returned if the estimated value goes
	  below \eqn{1-alpha}.}
	   \item{max.threads}{Decides the number of threads the program can use. Set to 0 for using the maximum number of threads allowed by the system (default).}
}
\value{
  A list:
  \item{P }{Value of the integral.}
  \item{Pe }{Estimated error of the P estimate.}
}
\author{
  David Bolin \email{davidbolin@gmail.com} and
  Finn Lindgren \email{finn.lindgren@gmail.com}
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
prob = excursions.int(mu.x, Q.x, mu.x-3, mu.x+3, max.threads=2)
prob$P
}
