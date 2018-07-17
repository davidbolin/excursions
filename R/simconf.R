## simconf.R
##
##   Copyright (C) 2012, 2013,2014 David Bolin, Finn Lindgren
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

#' Simultaneous confidence regions for Gaussian models
#'
#' \code{simconf} is used for calculating simultaneous confidence regions for
#' Gaussian models \eqn{x}. The function returns upper and lower bounds \eqn{a}
#' and \eqn{b} such that \eqn{P(a<x<b) = 1-alpha}.
#'
#' @param alpha Error probability for the region.
#' @param mu Expectation vector for the Gaussian distribution.
#' @param Q Precision matrix for the Gaussian distribution.
#' @param n.iter Number or iterations in the MC sampler that is used for
#' approximating probabilities. The default value is 10000.
#' @param Q.chol The Cholesky factor of the precision matrix (optional).
#' @param vars Precomputed marginal variances (optional).
#' @param ind Indices of the nodes that should be analyzed (optional).
#' @param verbose Set to TRUE for verbose mode (optional).
#' @param max.threads Decides the number of threads the program can use.
#' Set to 0 for using the maximum number of threads allowed by the system (default).
#' @param seed Random seed (optional).
#'
#' @return An object of class "excurobj" with elements
#' \item{a }{The lower bound.}
#' \item{b }{The upper bound.}
#' \item{a.marginal }{The lower bound for pointwise confidence bands.}
#' \item{b.marginal }{The upper bound for pointwise confidence bands.}
#' @export
#' @author David Bolin \email{davidbolin@gmail.com} and Finn Lindgren \email{finn.lindgren@gmail.com}
#' @references Bolin et al. (2015) \emph{Statistical prediction of global sea level
#' from global temperature}, Statistica Sinica, Vol 25, pp 351-367.
#' @seealso \code{\link{simconf.inla}}, \code{\link{simconf.mc}}, \code{\link{simconf.mixture}}
#' @examples
#' ## Create mean and a tridiagonal precision matrix
#' n = 11
#' mu.x = seq(-5, 5, length=n)
#' Q.x = Matrix(toeplitz(c(1, -0.1, rep(0, n-2))))
#' ## calculate the confidence region
#' conf = simconf(0.05, mu.x, Q.x, max.threads=2)
#' ## Plot the region
#' plot(mu.x, type="l", ylim=c(-10, 10),
#'      main='Mean (black) and confidence region (red)')
#' lines(conf$a, col=2)
#' lines(conf$b, col=2)

simconf <- function(alpha,
                    mu,
                    Q,
                    n.iter=10000,
                    Q.chol,
                    vars,
                    ind=NULL,
                    verbose=0,
                    max.threads=0,
                    seed=NULL)
{

  if(missing(mu)){
	  stop('Must specify mean value')
  } else {
    mu <- private.as.vector(mu)
  }
  if(missing(Q) && missing(Q.chol))
	  stop('Must specify a precision matrix or its Cholesky factor')

  if(!missing(ind) && !missing(Q.chol))
	  stop('Cannot provide both cholesky factor and indices.')

  if(!missing(Q))
    Q <- private.as.dgCMatrix(Q)

  if(!missing(Q.chol))
    Q.chol <- private.as.dtCMatrixU(Q.chol)

  if(!missing(vars))
    vars <- private.as.vector(vars)

  if(!missing(ind))
    ind <- private.as.vector(ind)


  if (!missing(Q.chol) && !is.null(Q.chol)) {
    L <- private.as.dtCMatrixU(Q.chol)
  } else {
    L <- suppressWarnings(private.Cholesky(Q, perm=FALSE)$R)
  }

  if(missing(vars)){
    vars  <- excursions.variances(L)
	}
	sd <- sqrt(vars)


  #setup function for optmization
  f.opt <- function(x,alpha,sd,L,ind,seed,max.threads){
	  q = qnorm(x)*sd;
	  prob = gaussint(a=-q, b=q, Q.chol=L, ind=ind, lim=1-alpha,
	 					  max.threads=max.threads,seed=seed)

	  if(prob$P == 0){
		  return(1)
	  } else {
		  return(prob$P)
	  }
  }

  r.o = optimize(f.opt,interval = c(0,1),alpha=alpha,sd=sd,L=L,seed=seed,max.threads = max.threads,ind=ind)

  a = mu-qnorm(r.o$minimum)*sd
  b = mu+qnorm(r.o$minimum)*sd

  a.marg = mu-qnorm(alpha/2)*sd
  b.marg = mu+qnorm(alpha/2)*sd


  if(is.null(ind)) {
    output <- list(a=a,b=b,a.marginal = a.marg,b.marginal=b.marg,
                   mean = mu, vars = vars)
  } else {
    output <- list(a=a[ind],b=b[ind],
                   a.marginal = a.marg[ind],
                   b.marginal=b.marg[ind],
                   mean = mu[ind], vars = vars[ind])
  }

  output$meta = list(calculation="simconf",
                     alpha=alpha,
                     n.iter=n.iter,
                     ind=ind,
                     call = match.call())
  class(output) <- "excurobj"
  return(output)
}

