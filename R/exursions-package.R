#' Excursions: Excursion Sets and Contour Credibility Regions for Random Fields
#'
#' \code{excursions} contains functions that compute probabilistic excursion sets,
#' contour credibility regions, contour avoiding regions, contour map quality measures,
#' and simultaneous confidence bands for latent Gaussian
#' random processes and fields. A detailed manual can be found in the paper 
#' Bolin, D and Lindgren, F  (2018) 
#' \emph{Calculating Probabilistic Excursion Sets and Related Quantities Using excursions},
#' Journal of Statistical Software, 86(5), 1--20.
#' 
#' The main functions in the package fall into three different categories described below.
#' 
#' \strong{Excursion sets, contour credibility regions, and contour avoiding regions}
#'
#' The main functions for computing excursion sets, contour credibility regions, and
#' contour avoiding regions are
#' \describe{
#' \item{\code{\link{excursions}} }{The main function for Gaussian models.}
#' \item{\code{\link{excursions.inla}} }{Interface for latent Gaussian models estimated using INLA.}
#' \item{\code{\link{excursions.mc}} }{Function for analyzing models that have been
#' estimated using Monte Carlo methods.}
#' }
#' The output from the functions above provides a discrete domain estimate of the regions.
#' Based on this estimate, the function \code{\link{continuous}} computes a continuous
#' domain estimate.
#'
#'  The main reference for these functions is Bolin, D. and Lindgren, F. (2015)
#'  \emph{Excursion and contour uncertainty regions for latent Gaussian models},
#'  JRSS-series B, vol 77, no 1, pp 85-106.
#'
#' \strong{Contour map quality measures}
#'
#' The package provides several functions for computing contour maps and their quality
#' measures. These quality measures can be used to decide on an appropriate number of
#' contours to use for the contour map.
#'
#' The main functions for computing contour maps and the corresponding quality measures
#' are
#' \describe{
#' \item{\code{\link{contourmap}} }{The main function for Gaussian models.}
#' \item{\code{\link{contourmap.inla}} }{Interface for latent Gaussian models estimated
#' using INLA.}
#' \item{\code{\link{contourmap.mc}} }{Function for analyzing models that have been
#' estimated using Monte Carlo methods.}
#' }
#' Other noteworthy functions relating to contourmaps are \code{\link{tricontour}} and
#' \code{\link{tricontourmap}}, which compute contour curves for functinos defined on
#' triangulations, as well as \code{\link{contourmap.colors}} which can be used to
#' compute appropriate colors for displaying contour maps.
#'
#' The main reference for these functions is Bolin, D. and Lindgren, F. (2017)
#' \emph{Quantifying the uncertainty of contour maps}, Journal of Computational and
#' Graphical Statistics, 26:3, 513-524.
#'
#' \strong{Simultaneous confidence bands}
#'
#' The main functions for computing simultaneous confidence bands are
#' \describe{
#' \item{\code{\link{simconf}} }{Function for analyzing Gaussian models.}
#' \item{\code{\link{simconf.inla}} }{Function for analyzing latent Gaussian models
#' estimated using INLA.}
#' \item{\code{\link{simconf.mc}} }{Function for analyzing models estimated using Monte
#' Carlo methods.}
#' \item{\code{\link{simconf.mixture}} }{Function for analyzing Gaussian mixture models.}
#' }
#'
#' The main reference for these functions is Bolin et al. (2015)
#' \emph{Statistical prediction of global sea level
#' from global temperature}, Statistica Sinica, Vol 25, pp 351-367.
#'
#' @importFrom graphics lines
#' @importFrom methods as is
#' @importFrom stats optimize pnorm qnorm quantile rnorm uniroot
#' @import Matrix
#' @import sp
#' @useDynLib excursions, .registration = TRUE
#' @aliases excursions-package
"_PACKAGE"