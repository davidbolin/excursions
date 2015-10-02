\name{contourmap.mc}
\alias{contourmap.mc}
\title{
  Contour maps and contour map quality measures using Monte Carlo samples
}
\description{
  \code{contourmap.mc} is used for calculating contour maps and quality measures for contour maps based on Monte Carlo samples of a model.
}
\usage{
contourmap.mc(samples,
              n.levels,
              ind,
              levels,
              type = c("standard",
                       "equalarea",
                       "P0-optimal",
                       "P1-optimal",
                       "P2-optimal"),
              compute = list(F=TRUE, measures = NULL),
              alpha,
              verbose=FALSE)
}
\arguments{
  \item{samples}{Matrix with model Monte Carlo samples. Each column contains a sample of the model.}
  \item{n.levels}{Number of levels in contour map}
  \item{ind}{Indices of the nodes that should be analyzed (optional)}
  \item{levels}{Levels to use in contour map}
  \item{type}{Type of contour map. One of:
    \itemize{
    \item{'standard' }{Equidistant levels between min(mu) and max(mu) (default)}
    \item{'equalarea' }{Levels such that different spatial regions are approximately equal in size}
    \item{'P0-optimal' }{Levels chosen to maximize the P0 measure}
    \item{'P1-optimal' }{Levels chosen to maximize the P1 measure}
    \item{'P2-optimal' }{Levels chosen to maximize the P2 measure}
    }}
  \item{compute}{A list with quality indices to compute
    \itemize{
      \item{'F': }{Compute the contour map function? (default TRUE)}
      \item{'measures': }{A list with the quality measures to compute
	("P0", "P1", or "P2")}
    }}
  \item{alpha}{Maximal error probability in contour map function (default=0.1)}
  \item{verbose}{Set to TRUE for verbose mode (optional)}
}
\value{
\code{contourmap} returns an object of class "excurobj". This is a list that can contains the following arguments:
  \item{u }{Contour levels}
  \item{n.levels}{The number of contours used}
  \item{u.e}{The values associated with the level sets G_k}
  \item{G}{A vector which shows which of the level sets G_k each node belongs to.}
  \item{map}{Representation of the contour map with map[i]=u.e[k] if i is in G_k}
  \item{F }{The contour map function (if computed)}
  \item{M }{Contour avoiding sets (if \code{F} is computed). \eqn{M=-1} for all non-significant nodes and  \eqn{M=k} for nodes that belong to \eqn{M_k}.}
  \item{P0/P1/P2}{Calculated quality measures (if computed)}
  \item{P0bound/P1bound/P2bound}{Calculated upper bounds quality measures (if computed)}
  \item{meta }{A list containing various information about the calculation.}
}
\author{
  David Bolin \email{davidbolin@gmail.com}
}

\references{
Bolin, D. and Lindgren, F. (2015) \emph{Quantifying the uncertainty of contour maps}, ArXiv preprint
}

\examples{
  n = 100
  Q = Matrix(toeplitz(c(1, -0.5, rep(0, n-2))))
  mu = seq(-5, 5, length=n)
  ## Sample the model 100 times (increase for better estimate)
  X = mu + solve(chol(Q),matrix(rnorm(n=n*100),nrow=n,ncol=100))

  lp <- contourmap.mc(X,n.levels = 2, compute=list(F=FALSE, measures = c("P1","P2")))
}
