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


gaussint <- function(mu,Q.chol, Q, a, b, lim = 0, n.iter = 10000, 
                     ind, max.size, max.threads=0,seed)
{

  if(missing(Q) && missing(Q.chol))
	  stop('Must specify a precision matrix or its Cholesky factor')

  if(missing(a))
    stop('Must specify lower integration limit')
    
  if(missing(b))
    stop('Must specify upper integration limit')
  
  if(!missing(ind) && !is.null(ind)){
    a[!ind] = Inf
    b[!ind] = -Inf
  }
  
  if (!missing(Q.chol) && !is.null(Q.chol)) {
      L = Q.chol
  } else {
      L = chol(Q)
  }

  if(!missing(mu) && !is.null(mu)){
    a = a - mu
    b = b - mu
  }   

  if (!is(L, "dtCMatrix"))
      stop("L needs to be in ccs format for now.")
    
  n = dim(L)[1]
  Mp = L@p
  Mi = L@i
  Mv = L@x
  
  if(missing(max.size))
    max.size = n
  

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

  return(list(Pv = out$Pv, Ev = out$Ev, P = out$Pv[1], E = out$Ev[1])) 
}
