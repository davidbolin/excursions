## excursions.R
##
##   Copyright (C) 2012, 2013, 2014, David Bolin, Finn Lindgren
##
##   This program is free software: you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation, either version 3 of the License, or
##   (at your option) any later version.
##
##   This program is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##
##   You should have received a copy of the GNU General Public License
##   along with this program.  If not, see <http://www.gnu.org/licenses/>.



# How to document packages:
# https://cran.r-project.org/web/packages/roxygen2/vignettes/rd.html
# In particular "this also works if there's already a function called pkgname()"
# However, that still leads to a name clash between excursions-package and
# excursions.  The @section and @rdname approach puts the documentation in
# the same file, with a separate heading for the package documentation.

#' @section Package:
#'
#' \code{excursions} contains functions that compute probabilistic excursion sets,
#' contour credibility regions, contour avoiding regions, contour map quality measures,
#' and simultaneous confidence bands for latent Gaussian
#' random processes and fields.
#' 
#' \strong{excursion sets, contour credibility regions, and contour avoiding regions}
#' 
#' The main functions for computing excursion sets, contour credibility regions, and 
#' contour avoiding regions are 
#' \itemize{
#' \item{\code{\link{excursions}} }{The main function for Gaussian models.}
#' \item{\code{\link{excursions.inla}} }{Interface for latent Gaussian models estimated using INLA.}
#' \item{\code{\link{excursions.mc}} }{Function for analyzing models that have been
#' estimated using Monte Carlo methods.}
#' }
#' The output from the functions above provides a discrete domain estimate of the regions. 
#' Based on this estimate, the function \code{\link{continuous}} computes a continuous 
#' domain estimate. 
#'  
#' \strong{Contour map quality measures}
#' 
#' The package provides several functinos for computing contour maps and their quality 
#' measures. These quality measures can be used to decide on an appropriate number of 
#' contours to use for the contour map.
#' 
#' The main functions for computing contour maps and the corresponding quality measures
#' are 
#' \itemize{
#' \item{\code{\link{contourmap}} }{The main function for Gaussian models.}
#' \item{\code{\link{contourmap.inla}} }{Interface for latent Gaussian models estimated using INLA.}
#' \item{\code{\link{contourmap.mc}} }{Function for analyzing models that have been
#' estimated using Monte Carlo methods.}
#' }
#' Other noteworthy functions relating to contourmaps are \code{\link{tricontour}} and
#' \code{\link{tricontourmap}}, which compute contour curves for functinos defined on
#' triangulations, as well as \code{\link{contourmap.colors}} which can be used to 
#' compute appropriate colors for displaying contour maps. 
#' 
#' \strong{Simultaneous confidence bands}
#' 
#' The main functions for computing simultaneous confidence bands are 
#' \itemize{\item{\code{\link{simconf}} }{Function for analyzing Gaussian models.}
#' \item{\code{\link{simconf.inla}} }{Function for analyzing latent Gaussian models
#' estimated using INLA.}
#' \item{\code{\link{simconf.mc}} }{Function for analyzing models estimated using Monte
#' Carlo methods.}
#' \item{\code{\link{simconf.mixture}} }{Function for analyzing Gaussian mixture models.}
#' }
#' 
#' 
#' @importFrom graphics lines
#' @importFrom methods as is
#' @importFrom stats optimize pnorm qnorm quantile rnorm uniroot
#' @import Matrix
#' @import sp
#' @useDynLib excursions, .registration = TRUE
#' 
#' @name excursions-package
#' @rdname excursions
#' 
NULL


#' Excursion Sets and Contour Credibility Regions for Random Fields
#' 
#' \code{excursions} is one of the main functions in the package with the same name. 
#' The function is used for calculating excursion sets, contour credible regions,
#' and contour avoiding sets for latent Gaussian models. Details on the function and the 
#' package are given in the sections below.
#'
#' @param alpha Error probability for the excursion set.
#' @param u Excursion or contour level.
#' @param mu Expectation vector.
#' @param Q Precision matrix.
#' @param type Type of region:
#'  \itemize{
#'     \item{'>' }{positive excursion region}
#'     \item{'<' }{negative excursion region}
#'     \item{'!=' }{contour avoiding region}
#'     \item{'=' }{contour credibility region}}
#' @param n.iter Number or iterations in the MC sampler that is used for approximating probabilities. The default value is 10000.
#' @param Q.chol The Cholesky factor of the precision matrix (optional).
#' @param F.limit The limit value for the computation of the F function. F is set to NA for all nodes where F<1-F.limit. Default is F.limit = \code{alpha}.
#' @param vars Precomputed marginal variances (optional).
#' @param rho Marginal excursion probabilities (optional). For contour regions, provide \eqn{P(X>u)}.
#' @param reo Reordering (optional).
#' @param method Method for handeling the latent Gaussian structure:
#'  \itemize{
#'       \item{'EB' }{Empirical Bayes (default)}
#'       \item{'QC' }{Quantile correction, rho must be provided if QC is used.}}
#' @param ind Indices of the nodes that should be analysed (optional).
#' @param max.size Maximum number of nodes to include in the set of interest (optional).
#' @param verbose Set to TRUE for verbose mode (optional).
#' @param max.threads Decides the number of threads the program can use. Set to 0 for using the maximum number of threads allowed by the system (default).
#' @param seed Random seed (optional).
#'
#' @return \code{excursions} returns an object of class "excurobj". This is a list that contains the following arguments:
#' \item{E }{Excursion set, contour credible region, or contour avoiding set}
#' \item{G }{ Contour map set. \eqn{G=1} for all nodes where the \eqn{mu > u}.}
#' \item{M }{ Contour avoiding set. \eqn{M=-1} for all non-significant nodes. \eqn{M=0} for nodes where the process is significantly below \code{u} and \eqn{M=1} for all nodes where the field is significantly above \code{u}. Which values that should be present depends on what type of set that is calculated.}
#' \item{F }{The excursion function corresponding to the set \code{E} calculated or values up to \code{F.limit}}
#' \item{rho }{Marginal excursion probabilities}
#' \item{mean }{The mean \code{mu}.}
#' \item{vars }{Marginal variances.}
#' \item{meta }{A list containing various information about the calculation.}
#' @export
#' @details 
#' The estimation of the region is done using sequential importance sampling with 
#' \code{n.iter} samples. The procedure requires computing the marginal variances of 
#' the field, which should be supplied if available. If not, they are computed using 
#' the Cholesky factor of the precision matrix. The cost of this step can therefore be 
#' reduced by supplying the Cholesky factor if it is available.
#' 
#' The latent structure in the latent Gaussian model can be handled in several different 
#' ways. The default strategy is the EB method, which is
#' exact for problems with Gaussian posterior distributions. For problems with
#' non-Gaussian posteriors, the QC method can be used for improved results. In order to use 
#' the QC method, the true marginal excursion probabilities must be supplied using the 
#' argument \code{rho}. 
#' Other more
#' complicated methods for handling non-Gaussian posteriors must be implemented manually
#' unless \code{INLA} is used to fit the model. If the model is fitted using \code{INLA},
#' the method \code{excursions.inla} can be used. See the Package section for further details
#' about the different options.
#' @author David Bolin \email{davidbolin@@gmail.com} and Finn Lindgren \email{finn.lindgren@@gmail.com}
#' @references Bolin, D. and Lindgren, F. (2015) \emph{Excursion and contour uncertainty regions for latent Gaussian models}, JRSS-series B, vol 77, no 1, pp 85-106.
#' @seealso \code{\link{excursions.inla}}, \code{\link{excursions.mc}}
#'
#' @examples
#' ## Create a tridiagonal precision matrix
#' n = 21
#' Q.x = sparseMatrix(i=c(1:n, 2:n), j=c(1:n, 1:(n-1)), x=c(rep(1, n), rep(-0.1, n-1)),
#'                    dims=c(n, n), symmetric=TRUE)
#' ## Set the mean value function
#' mu.x = seq(-5, 5, length=n)
#' 
#' ## calculate the level 0 positive excursion function
#' res.x = excursions(alpha=1, u=0, mu=mu.x, Q=Q.x, 
#'                    type='>', verbose=1, max.threads=2)
#'                    
#' ## Plot the excursion function and the marginal excursion probabilities
#' plot(res.x$F, type="l", 
#'      main='Excursion function (black) and marginal probabilites (red)')
#' lines(res.x$rho, col=2)

excursions <- function(alpha,
                       u,
                       mu,
                       Q,
                       type,
                       n.iter=10000,
                       Q.chol,
                       F.limit,
                       vars,
                       rho,
                       reo,
                       method='EB',
                       ind,
                       max.size,
                       verbose=0,
                       max.threads=0,
                       seed)
{

  if(method=='QC'){
    qc = TRUE
  } else if(method == 'EB'){
    qc = FALSE
  } else {
    stop('only EB and QC methods are supported.')
  }
  if(missing(alpha))
    stop('Must specify error probability')

  if(missing(u))
    stop('Must specify level')

  if(missing(mu)){
    stop('Must specify mean value')
  } else {
    mu <- private.as.vector(mu)
  }
  if(missing(Q) && missing(Q.chol))
    stop('Must specify a precision matrix or its Cholesky factor')

  if(missing(type))
    stop('Must specify type of excursion set')

  if(qc && missing(rho))
    stop('rho must be provided if QC is used.')

  if(!missing(ind) && !missing(reo))
    stop('Either provide a reordering using the reo argument or provied a set of nodes using the ind argument, both cannot be provided')

  if(missing(F.limit)) {
    F.limit = alpha
  } else {
    F.limit = max(alpha,F.limit)
  }

  if (!missing(Q.chol) && !is.null(Q.chol)) {
    ## make the representation unique (i,j,v) and upper triangular
    Q = private.as.dgTMatrix(private.as.dtCMatrixU(Q.chol))
    is.chol = TRUE
  } else {
    ## make the representation unique (i,j,v)
    Q = private.as.dgTMatrix(Q)
    is.chol = FALSE
  }

  if (missing(vars)) {
    if(is.chol){
      vars <- excursions.variances(L=Q)
    } else {
      vars <- excursions.variances(Q=Q)
    }
  } else {
    vars <- private.as.vector(vars)
  }

  if(!missing(rho))
    rho <- private.as.vector(rho)

  if(!missing(ind))
    ind <- private.as.vector(ind)


  if(verbose)
    cat("Calculate marginals\n")
  marg <- excursions.marginals(type = type, rho = rho,vars = vars,
                               mu = mu, u = u, QC = qc)

  if (missing(max.size)){
    m.size = length(mu)
  } else {
    m.size = max.size
  }
  if (!missing(ind)) {
    if(is.logical(ind)){
      indices = ind
      if(missing(max.size)){
        m.size = sum(ind)
      } else {
        m.size = min(sum(ind),m.size)
      }
    } else {
      indices = rep(FALSE,length(mu))
      indices[ind] = TRUE
      if(missing(max.size)){
        m.size = length(ind)
      } else {
        m.size = min(length(ind),m.size)
      }
    }
  } else {
    indices = rep(TRUE,length(mu))
  }

  if(verbose)
    cat("Calculate permutation\n")
  if(missing(reo)){
    use.camd = !missing(ind) || F.limit < 1
    if(qc){
      reo <- excursions.permutation(marg$rho_ng, indices,
                                    use.camd = TRUE,F.limit,Q)
    } else {
      reo <- excursions.permutation(marg$rho, indices,
                                    use.camd = TRUE,F.limit,Q)
    }
  } else {
    reo <- private.as.vector(reo)
  }

  if(verbose)
    cat("Calculate limits\n")
  limits <- excursions.setlimits(marg, vars,type,QC=qc,u,mu)

  res <- excursions.call(limits$a,limits$b,reo,Q, is.chol = is.chol,
                         1-F.limit, K = n.iter, max.size = m.size,
                         n.threads = max.threads,seed = seed)

  n = length(mu)
  ii = which(res$Pv[1:n] > 0)
  if (length(ii) == 0) i=n+1 else i=min(ii)

  F = Fe  = E = G = rep(0,n)
  F[reo] = res$Pv
  Fe[reo] = res$Ev

  ireo = NULL
  ireo[reo] = 1:n

  ind.lowF = F < 1-F.limit
  E[F>1-alpha] = 1

  if(type == '=') {
    F=1-F
  }

  if(type == "<") {
    G[mu>u] = 1
  } else {
    G[mu>=u] = 1
  }

  F[ind.lowF] = Fe[ind.lowF] = NA

  M = rep(-1,n)
  if (type=="<") {
    M[E==1] = 0
  } else if (type == ">") {
    M[E==1] = 1
  } else if (type == "!=" || type == "=") {
    M[E==1 & mu>u] = 1
    M[E==1 & mu<u] = 0
  }

  if (missing(ind) || is.null(ind)) {
    ind <- seq_len(n)
  } else if(is.logical(ind)){
    ind <- which(ind)
  }

  output <- list(F = F,
                 G = G,
                 M = M,
                 E = E,
                 mean = mu,
                 vars=vars,
                 rho=marg$rho,
                 meta=(list(calculation="excursions",
                            type=type,
                            level=u,
                            F.limit=F.limit,
                            alpha=alpha,
                            n.iter=n.iter,
                            method=method,
                            ind=ind,
                            reo=reo,
                            ireo=ireo,
                            Fe=Fe,
                            call = match.call())))
  class(output) <- "excurobj"
  output
}
