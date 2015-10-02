\name{contourmap.colors}
\alias{contourmap.colors}
\title{
  Define a color map for displaying contour maps.
}
\description{
  \code{contourmap.colors} calculates suitable colours for displaying contour maps.
}
\usage{
contourmap.colors(lp,
                  zlim,
                  col,
                  credible.col)
}
\arguments{
  \item{lp }{A contourmap calculated by contourmap.}
  \item{zlim }{The range that should be used (optional). The default is the range of the mean value function used when creating the contourmap.}
  \item{col }{The colormap that the colours should be taken from.}
  \item{credible.col }{The color that should be used for displaying the credible regions for the contour curves (optional).}
}
\value{
  A color map.}
\author{
  David Bolin \email{davidbolin@gmail.com}
}

\examples{
  n = 10
  Q = Matrix(toeplitz(c(1, -0.5, rep(0, n-2))))
  map <- contourmap(mu = seq(-5, 5, length=n),Q,n.levels = 2,compute=list(F=FALSE),max.threads=1)
  cols = contourmap.colors(map, col=heat.colors(100, 1),
                           credible.col = grey(0.5, 1))
}
