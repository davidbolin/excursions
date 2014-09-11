\name{contourmap}
\alias{contourmap}
\title{
  Contour maps and contour map quality measures for latent Gaussian models
}
\description{
  \code{contourmap} is used for calculating contour maps and quality measures for contour maps for latent Gaussian models.
}
\usage{
contourmap <- function(mu,
                       Q,
                       vars,
                       n.levels,
                       ind,
                       levels,
                       type='standard',
                       measure,
                       contour.map.function = FALSE,
                       use.marginals=TRUE,
                       alpha=1,
                       n.iter=10000,
                       verbose=FALSE,
                       max.threads=0)
}
\arguments{
  \item{mu}{Expectation vector}
  \item{Q}{Precision matrix}
  \item{vars}{Precomputed marginal variances (optional)}
  \item{n.levels}{Number of levels in contour map}
  \item{ind}{Indices of the nodes that should be analyzed (optional)}
  \item{levels}{Levels to use in contour map}
  \item{type}{Type of contour map:
    \itemize{
    \item{'standard' }{Equidistant levels between min(mu) and max(mu) (default)}
    \item{'equalarea' }{Levels such that different spatial regions are approximately equal in size}
    \item{'P0-optimal' }{Levels chosen to maximize the P0 measure}
    \item{'P1-optimal' }{Levels chosen to maximize the P1 measure}
    \item{'P2-optimal' }{Levels chosen to maximize the P2 measure}
    }}
  \item{measure}{Vector of quality measures to compute, measures 0, 1, and 2 are available (optional).}
  \item{contour.map.function}{Indicates whether the contour map function should be calculated (default=FALSE)}
  \item{alpha}{Maximal error probability in contour map function (default=1)}
  \item{n.iter}{Number or iterations in the MC sampler that is used for
    approximating probabilities. The default value is 10000.}
  \item{verbose}{Set to TRUE for verbose mode (optional)}
  \item{max.threads}{Decides the number of threads the program can use. Set to 0 for using the maximum number of threads allowed by the system (default).}
}
\value{
  A list:
  \item{u }{Contour levels}
  \item{n.levels}{The number of contours used}
  \item{u.e}{The values associated with the level sets E_k}
  \item{E}{A list containing the different level sets E_k}
  \item{map}{Representation of the contour map with map[i]=u.e[k] if i is in E[[k]]}
  \item{F }{The contour map function}
  \item{P0/P1/P2}{Calculated quality measures}
  \item{P0bound/P1bound/P2bound}{Calculated upper bounds quality measures}
}
\author{
  David Bolin \email{davidbolin@gmail.com}
}

\references{
Bolin, D. and Lindgren, F. (2014) \emph{Quantifying the uncertainty of contour maps}, ArXiv preprint
}

\examples{
  n = 100
  Q = Matrix(toeplitz(c(1, -0.5, rep(0, n-2))))
  mu = seq(-5, 5, length=n)
  lp <- contourmap(mu,Q,n.levels = 2, measure = c(1,2))
}
