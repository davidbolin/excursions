## excursions.int.R
##
##   Copyright (C) 2012, 2013, 2014 David Bolin, Finn Lindgren
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

#' Sequential estimation of Gaussian integrals
#'
#' \code{gaussint} is used for calculating \eqn{n}-dimensional Gaussian integrals
#' \deqn{\int_a^b \frac{|Q|^{1/2}}{(2\pi)^{n/2}}
#' \exp(-\frac1{2}(x-\mu)^{T}Q(x-\mu)) dx}{|Q|^(1/2)*(2\pi)^(-n/2) \int_a^b exp(-0.5*(x-\mu)^T Q (x-\mu)) dx}
#' A limit value \eqn{lim} can be used to stop the integration if the sequential
#' estimate goes below the limit, which can result in substantial computational
#' savings in cases when one only is interested in testing if the integral is above
#' the limit value. The integral is calculated sequentially, and estimates for
#' all subintegrals are also returned.
#'
#' @param mu Expectation vector for the Gaussian distribution.
#' @param Q.chol The Cholesky factor of the precision matrix (optional).
#' @param Q Precision matrix for the Gaussian distribution. If Q is supplied but not Q.chol,
#' the cholesky factor is computed before integrating.
#' @param a Lower limit in integral.
#' @param b Upper limit in integral.
#' @param lim If this argument is used, the integration is stopped and 0 is returned
#' if the estimated value goes below \eqn{lim}.
#' @param n.iter Number or iterations in the MC sampler that is used for approximating
#' probabilities. The default value is 10000.
#' @param ind Indices of the nodes that should be analyzed (optional).
#' @param use.reordering Determines what reordering to use:
#' \itemize{
#'   \item{"natural" }{No reordering is performed.}
#'   \item{"sparsity" }{Reorder for sparsity in the cholesky factor (MMD reordering
#'   is used).}
#'   \item{"limits" }{Reorder by moving all nodes with a=-Inf and b=Inf first and
#'   then reordering for sparsity (CAMD reordering is used).}}
#' @param max.size The largest number of sub-integrals to compute. Default is the total
#' dimension of the distribution.
#' @param max.threads Decides the number of threads the program can use. Set to 0 for
#' using the maximum number of threads allowed by the system (default).
#' @param seed The random seed to use (optional).
#'
#' @return A list with elements
#' \item{P }{Value of the integral.}
#' \item{E }{Estimated error of the P estimate.}
#' \item{Pv }{A vector with the estimates of all sub-integrals.}
#' \item{Ev }{A vector with the estimated errors of the Pv estimates.}
#' @export
#' @author David Bolin \email{davidbolin@gmail.com}
#' @references Bolin, D. and Lindgren, F. (2015) \emph{Excursion and contour uncertainty regions for latent Gaussian models}, JRSS-series B, vol 77, no 1, pp 85-106.
#' @examples
#' ## Create mean and a tridiagonal precision matrix
#' n = 11
#' mu.x = seq(-5, 5, length=n)
#' Q.x = Matrix(toeplitz(c(1, -0.1, rep(0, n-2))))
#' ## Calculate the probability that the variable is between mu-3 and mu+3
#' prob = gaussint(mu=mu.x, Q=Q.x, a= mu.x-3, b=mu.x+3, max.threads=2)
#' prob$P

gaussint <- function(mu,
                     Q.chol,
                     Q,
                     a,
                     b,
                     lim = 0,
                     n.iter = 10000,
                     ind,
                     use.reordering = c("natural","sparsity","limits"),
                     max.size,
                     max.threads=0,
                     seed)
{

  if( missing(Q) && missing(Q.chol))
	  stop('Must specify a precision matrix or its Cholesky factor')

  if(missing(a))
    stop('Must specify lower integration limit')

  if(missing(b))
    stop('Must specify upper integration limit')

  if(!missing(mu))
    mu <- private.as.vector(mu)

  if(!missing(ind))
    ind <- private.as.vector(ind)

  a <- private.as.vector(a)
  b <- private.as.vector(b)

  if(!missing(Q))
    Q <- private.as.dgCMatrix(Q)

  if(!missing(Q.chol))
    Q.chol <- private.as.dtCMatrixU(Q.chol)


  n = length(a)
  if(length(b) != n)
    stop('Vectors with integration limits are of different length.')

  use.reordering <- match.arg(use.reordering)

  if(!missing(ind) && !is.null(ind)){
    a[!ind] = -Inf
    b[!ind] = Inf
  }

  if(missing(max.size))
    max.size = n

  reordered = FALSE
  if(!missing(Q.chol) && !is.null(Q.chol)){
    ## Cholesky factor is provided, use that and do not reorder
      L = Q.chol
      if(dim(L)[1] != dim(L)[2]){
        stop("Q.chol is not square")
      } else if(dim(L)[1] != n) {
        stop("Dimensions of Q.chol is different from the length of the integration limits.")
      }
  } else if(!missing(Q) && !is.null(Q)){
    if(dim(Q)[1] != dim(Q)[2]){
      stop("Q is not square")
      } else if(dim(Q)[1] != n) {
        stop("Dimensions of Q is different from the length of the integration limits.")
      }
    ## Cholesky factor is not provided and we are told to reorder
    if(use.reordering == "limits")
    {
      #Reorder by moving nodes with limits -inf ... inf first
      inf.ind = (a==-Inf) & (b==Inf)
      if(sum(inf.ind)>0){
        max.size = min(max.size,n-sum(inf.ind))
        n = length(a)
        cind = rep(1,n)
        cind[inf.ind] = 0
        reo = rep(0,n)
        out <- .C("reordering",nin = as.integer(n), Mp = as.integer(Q@p),
                  Mi = as.integer(Q@i), reo = as.integer(reo),
                  cind = as.integer(cind))
        reo = out$reo+1
        Q = Q[reo,reo]
        reordered = TRUE
      }

      L <- suppressWarnings(private.Cholesky(Q,perm=FALSE)$R)

    } else if(use.reordering == "sparsity"){
      ## Reorder for sparsity calculated by Matrix
      L <- suppressWarnings(private.Cholesky(Q,perm=TRUE))
      reo <- L$reo
      ireo <- L$ireo
      reordered = TRUE
      L <- L$R
    } else {
      ## Do not reorder
      L <- suppressWarnings(private.Cholesky(Q,perm=FALSE)$R)
    }
  }

  # Note: If lim > 0 and reorder == TRUE, we should calculate marginal
  # probabilities, see if bound is above lim, and then reorder

  if(!missing(mu) && !is.null(mu)){
    if(length(mu) != n){
      stop("The length of mu is different from the length of the integration limits.")
    }
    a = a - mu
    b = b - mu
  }

  a[a==Inf]  = .Machine$double.xmax
  b[b==Inf]  = .Machine$double.xmax
  a[a==-Inf] = -.Machine$double.xmax
  b[b==-Inf] = -.Machine$double.xmax

	if(reordered == TRUE){
		a = a[reo]
		b = b[reo]
  }

  if (is(L, "Matrix")) {
     L <- private.as.dtCMatrixU(L)
  } else {
    stop("Unsuported matrix type.")
  }

  if(!missing(seed) && !is.null(seed)){
    seed_provided = 1
    seed <- as.integer(seed)
    if(length(seed)==1){
      seed <- rep(seed,6)
    }
    seed.in = seed
  } else {
    seed_provided = 0
    seed.in = as.integer(rep(0,6))
  }

  Pv = Ev = rep(0,dim(L)[1])

  opts = c(n,n.iter,max.size,max.threads,seed_provided)

  out <- .C("shapeInt", Mp = as.integer(L@p), Mi = as.integer(L@i),
              Mv = as.double(L@x), a = as.double(a), b = as.double(b),
              opts = as.integer(opts), lim = as.double(lim),
              Pv = as.double(Pv), Ev = as.double(Ev),seed_in=seed.in)

  P = out$Pv[1]
  E = out$Ev[1]
  if(use.reordering == "limits") {
    out$Pv[1:(dim(L)[1]-max.size)] = out$Pv[dim(L)[1]-max.size+1]
    out$Ev[1:(dim(L)[1]-max.size)] = out$Ev[dim(L)[1]-max.size+1]
    P = out$Pv[1]
    E = out$Ev[1]
  } else if(use.reordering == "sparsity") {
    out$Pv = out$Pv[ireo]
    out$Ev = out$Ev[ireo]
  }
  return(list(Pv = out$Pv, Ev = out$Ev, P = P, E = E))
}
