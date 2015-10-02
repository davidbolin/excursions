\name{submesh.grid}
\alias{submesh.grid}
\title{
  Extract a part of a grid
}
\description{
  Extracts a part of a grid
}
\usage{
submesh.grid(z, grid=NULL)
}
\arguments{
  \item{z}{A matrix with values indicating which nodes that should be present in the submesh.}
\item{grid}{A list with locations and dimensions of the grid.}
}

\value{
  An \code{inla.mesh} object.
}
\note{
  This function requires the \code{INLA} package, which is not a CRAN
  package.  See \url{http://www.r-inla.org/download} for easy
  installation instructions.
}

\author{
Finn Lindgren
}
\examples{
\dontrun{
if (require(INLA)) {
  nxy = 40
  x=seq(from=0,to=4,length.out=nxy)
  lattice=inla.mesh.lattice(x=x,y=x)
  mesh=inla.mesh.create(lattice=lattice, extend=FALSE, refine=FALSE)

  #extract a part of the mesh inside a circle
  xy.in <- rowSums((mesh$loc[,1:2]-2)^2)<1
  submesh <- excursionsdevel:::submesh.grid(matrix(xy.in,nxy,nxy),
                                          list(loc=mesh$loc,dim=c(nxy,nxy)))
  plot(mesh$loc[,1:2])
  lines(2+cos(seq(0,2*pi,length.out=100)), 2+sin(seq(0,2*pi,length.out=100)))
  plot(submesh,add=TRUE)
  points(mesh$loc[xy.in,1:2],col="2")
}
}
}