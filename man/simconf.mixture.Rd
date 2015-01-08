\name{simconf.mixture}
\alias{simconf.mixture}
\title{
  Simultaneous confidence regions for Gaussian mixture models
}
\description{
  \code{excursions.simconf} is used for calculating simultaneous
  confidence regions for Gaussian mixture models. The distribution for the process \eqn{x} is assumed to be
 \deqn{latex}{
  \pi(x) = \sum_{k=1}^K w_k N(\mu_k, Q_k^{-1}).
 }
The function returns upper and lower bounds \eqn{a} and \eqn{b} such that \eqn{P(a<x<b) = 1-alpha}.
}
\usage{
simconf.mixture(alpha,
                mu,
                Q,
                w,
                ind,
                n.iter=10000,
                vars,
                verbose=0,
                max.threads=0,
                seed=NULL,
                mix.samp = TRUE)
}
\arguments{
  \item{alpha }{Error probability for the region.}
  \item{mu }{A list with the \code{k} expectation vectors \deqn{latex}{\mu_k}.}
  \item{Q }{A list with the \code{k} precision matrices \deqn{latex}{Q_k}.}
  \item{w }{A vector with the weights for each class in the mixture.}
  \item{ind}{Indices of the nodes that should be analyzed
    (optional)}
  \item{n.iter}{Number or iterations in the MC sampler that is
    used for approximating probabilities. The default value is
    10000.}
  \item{vars }{A list with precomputed marginal variances for each class (optional).}
  \item{verbose }{Set to TRUE for verbose mode (optional)}
  \item{max.threads }{Decides the number of threads the program can use. Set to 0 for using the maximum number of threads allowed by the system (default).}
  \item{seed }{Random seed (optional)}
  \item{mix.samp }{If TRUE, do MC integration by directly sampling the mixture, otherwise use sequential integration.}

}
\value{
  A list:
  \item{a }{The lower bound.}
  \item{b }{The upper bound.}
  \item{a.marginal }{The lower bound for pointwise confidence bands.}
  \item{b.marginal }{The upper bound for pointwise confidence bands.}

}
\author{
  David Bolin \email{davidbolin@gmail.com}
}

\references{
Bolin et al. Statistical prediction of global sea level from global temperature, Statistica Sinica, Vol 25, pp 351-367.
}

\examples{
n = 11
K = 3
mu <- Q <- list()
for(k in 1:K){
  mu[[k]] = k*0.1 + seq(-5, 5, length=n)
  Q[[k]] = Matrix(toeplitz(c(1, -0.1, rep(0, n-2))))
}

## calculate the confidence region
conf = simconf.mixture(0.05, mu, Q, w = rep(1/3,3), max.threads=2)

## Plot the region
plot(mu[[1]],type="l")
lines(mu[[2]])
lines(mu[[3]])
lines(conf$a, col=2)
lines(conf$b, col=2)
}
