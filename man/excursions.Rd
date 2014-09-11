\name{excursions}
\alias{excursions}
\title{
  Excursion sets and contour credible regions for latent Gaussian
models
}
\description{
  \code{excursions} is used for calculating excursion sets, contour
  credible regions, and contour avoiding sets for latent Gaussian
  models. The latent structure can be handled in several different
  ways. The default strategy is the EB method, which is exact for
  problems with Gaussian posterior distributions. For problems with
  non-Gaussian posteriors, the QC method can be used for improved
  results. Other more complicated methods for handling non-Gaussian
  posteriors must be implemented manually unless \code{INLA} is used to
  fit the model. If the model is fitted using \code{INLA}, the method
  \code{excursions.inla} can be used.
}
\usage{
excursions <- function(alpha,
                       u,
                       mu,
                       Q,
                       type,
                       n.iter=10000,
                       Q.chol,
                       vars,
                       rho,
                       reo,
                       method='EB',
                       ind,
                       max.size,
                       verbose=0,
                       max.threads=0,
                       seed)
}
\arguments{
  \item{alpha}{Error probability. The excursion function is calculated
    for all values up to \code{alpha}. }
  \item{u}{Excursion or contour level.}
  \item{mu}{Expectation vector}
  \item{Q}{Precision matrix}
  \item{type}{Type of region:
    \itemize{
    \item{'>' }{positive excursion region}
    \item{'<' }{negative excursion region}
    \item{'!=' }{contour avoiding region}
    \item{'=' }{contour credibility region}}}
  \item{n.iter}{Number or iterations in the MC sampler that is used for
    approximating probabilities. The default value is 10000.}
  \item{Q.chol}{The Cholesky factor of the precision matrix (optional)}
  \item{vars}{Precomputed marginal variances (optional)}
  \item{rho}{Marginal excursion probabilities (optional).  For contour
    regions, provide \eqn{P(X>u)}.}
  \item{reo}{Reordering (optional)}
  \item{ind}{Indices of the nodes that should be analyzed (optional)}
  \item{max.size}{Maximum number of nodes to include in the set of
    interest (optional)}
  \item{method}{Method for handeling the latent Gaussian structure:
  \itemize{
    \item{'EB' }{Empirical Bayes (default)}
    \item{'QC' }{Quantile correction, rho must be provided if QC is used.}
  }
  }
  \item{verbose}{Set to TRUE for verbose mode (optional)}
  \item{max.threads}{Decides the number of threads the program can use. Set to 0 for using the maximum number of threads allowed by the system (default).}
}
\value{
  A list:
  \item{D }{Excursion set, contour credible region, or contour avoiding
    set}
  \item{F }{The excursion function corresponding to the set D calculated
    for values up to alpha}
  \item{Fe }{Estimated standard deviations of F estimate}
  \item{rho }{Marginal excursion probabilities}
  \item{reo }{Reordering used for excursion function}
  \item{ireo }{Inverse reordering}
  \item{vars }{Marginal variances}
}
\author{
  David Bolin \email{davidbolin@gmail.com} and
  Finn Lindgren \email{finn.lindgren@gmail.com}
}

\references{
Bolin, D. and Lindgren, F. (2013) \emph{Excursion and contour uncertainty regions for latent Gaussian models}, Journal of the Royal Statistical Society: Series B (in press).
}

\examples{
## Create a tridiagonal precision matrix
n = 21
Q.x = sparseMatrix(i=c(1:n, 2:n),
                   j=c(1:n, 1:(n-1)),
                   x=c(rep(1, n), rep(-0.1, n-1)),
                   dims=c(n, n),
                   symmetric=TRUE)

## Set the mean value function
mu.x = seq(-5, 5, length=n)

## calculate the level 0 positive excursion function
res.x = excursions(alpha=1, u=0, mu=mu.x, Q=Q.x, type='>',
                   verbose=1, max.threads=2)

## Plot the excursion function and the marginal excursion probabilities
plot(res.x$F, type="l",
     main='Excursion function (black) and marginal probabilites (red)')
lines(res.x$rho, col=2)
}
