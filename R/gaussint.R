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


gaussint <- function(mu, Q.chol, Q, a, b, lim = 0, n.iter = 10000,
                     ind, max.size, max.threads=0,seed, reo = FALSE)
{

  if( missing(Q) && missing(Q.chol))
	  stop('Must specify a precision matrix or its Cholesky factor')

  if(missing(a))
    stop('Must specify lower integration limit')

  if(missing(b))
    stop('Must specify upper integration limit')

  if(!missing(ind) && !is.null(ind)){
    a[!ind] = Inf
    b[!ind] = -Inf
  }
  n = length(a)
  if(missing(max.size))
    max.size = n

  reordered = FALSE
  if(reo && !missing(Q) && !is.null(Q)){
    inf.ind = (a==-Inf) & (b==Inf)
    if(sum(inf.ind)>0){
      # We have nodes with limits -inf .. inf
      # move these first in vector
      max.size = min(max.size,n-sum(inf.ind))
      n = length(a)
      cind = rep(1,n)
      cind[inf.ind] = 0
      reo = rep(0,n)
			out <- .C("reordering",nin = as.integer(n), Mp = as.integer(Q@p),
		                        Mi = as.integer(Q@i), reo = as.integer(reo),
		                        cind = as.integer(cind))
		  reo = out$reo+1
		  a = a[reo]
		  b = b[reo]
		  Q = Q[reo,reo]
		  reordered = TRUE
    }
    L = chol(Q)
  } else {
    if(!missing(Q.chol) && !is.null(Q.chol))
      L = Q.chol
    else
     L = chol(Q)
  }
    # If lim > 0 and reorder == TRUE, we should calculate marginal
    # probabilities, see if bound is above lim, and then eorder


  if(!missing(mu) && !is.null(mu)){
    a = a - mu
    b = b - mu
  }

  a[a==Inf]  = .Machine$double.xmax
  b[b==Inf]  = .Machine$double.xmax
  a[a==-Inf] = -.Machine$double.xmax
  b[b==-Inf] = -.Machine$double.xmax

  if (!is(L, "dtCMatrix"))
      stop("L needs to be in ccs format for now.")

  n = dim(L)[1]
  Mp = L@p
  Mi = L@i
  Mv = L@x

  if(!missing(seed) && !is.null(seed)){
    seed_provided = 1
    seed.in = seed
  } else {
    seed_provided = 0
    seed.in = as.integer(rep(0,6))
  }

  Pv = rep(0,n)
  Ev = rep(0,n)

  opts = c(n,n.iter,max.size,max.threads,seed_provided)

  out <- .C("shapeInt", Mp = as.integer(Mp), Mi = as.integer(Mi),
              Mv = as.double(Mv), a = as.double(a), b = as.double(b),
              opts = as.integer(opts), lim = as.double(lim),
              Pv = as.double(Pv), Ev = as.double(Ev),seed_in=seed.in)

  if(reordered) {
    out$Pv[1:(n-max.size)] = out$Pv[n-max.size+1]
    out$Ev[1:(n-max.size)] = out$Ev[n-max.size+1]
  }
  return(list(Pv = out$Pv, Ev = out$Ev, P = out$Pv[1], E = out$Ev[1]))
}
