\name{simconf.inla}
\alias{simconf.inla}
\title{
  Simultaneous confidence regions for latent Gaussian models
}
\description{
  \code{simconf.inla} is used for calculating simultaneous
  confidence regions for latent Gaussian models estimated using \code{INLA}.
}
\usage{
simconf.inla(result.inla,
             stack,
             name=NULL,
             tag=NULL,
             ind=NULL,
             alpha,
             method="NI",
             n.iter=10000,
             verbose=0,
             max.threads=0,
             seed=NULL,
             inla.sample=TRUE)
}
\arguments{
  \item{result.inla }{Result object from INLA call}
  \item{stack }{The stack object used in the INLA call.}
  \item{tag }{The tag of the component in the stack for which to do the calculation. This argument should only be used if a stack object is provided, use the name argument otherwise.}
  \item{name }{The name of the component for which to do the calculation. This argument should only be used if a stack object is not provided, use the tag argument otherwise.}
  \item{ind }{If only a part of a component should be used in the calculations, this argument specifies the indices for that part.}
  \item{method }{Method for handeling the latent Gaussian structure:
    \itemize{
      \item{'EB' }{Empirical Bayes (Gaussian approximation of posterior)}
      \item{'NI' }{Numerical integration (Calculation based on the Gaussian mixture approximation of the posterior, as calculated by INLA)}
    }
  }
  \item{alpha }{Error probability for the region.}
  \item{n.iter }{Number or iterations in the MC sampler that is
    used for approximating probabilities. The default value is
    10000.}
  \item{verbose }{Set to TRUE for verbose mode (optional)}
  \item{max.threads }{Decides the number of threads the program can use. Set to 0 for using the maximum number of threads allowed by the system (default).}
  \item{seed }{Random seed (optional)}
  \item{inla.sample }{Use inla.posterior.sample for MC integration?}
}
\value{
  A list:
  \item{a }{The lower bound.}
  \item{b }{The upper bound.}
  \item{a.marginal }{The lower bound for pointwise confidence bands.}
  \item{b.marginal }{The upper bound for pointwise confidence bands.}

}
\note{
  This function requires the \code{INLA} package, which is not a CRAN
  package.  See \url{http://www.r-inla.org/download} for easy
  installation instructions.
}
\author{
  David Bolin \email{davidbolin@gmail.com}
}

\examples{
if (require(INLA)) {
n <- 10
x <- seq(0, 6, length.out=n)
y <- sin(x) + rnorm(n)
mu <- 1:n

result <- inla(y ~ 1 + f(mu, model='rw2'),
               data=list(y=y, mu=mu), verbose=FALSE,
               control.compute = list(config=TRUE))

res <- simconf.inla(result, name='mu', alpha = 0.05, max.threads=2)

plot(result$summary.random$mu$mean,ylim=c(-2,2))
lines(res$a)
lines(res$b)
lines(res$a.marginal,col="2")
lines(res$b.marginal,col="2")
}
}