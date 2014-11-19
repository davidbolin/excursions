## excursions.simconf.R
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

simconf <- function(alpha, mu, Q, n.iter=10000, Q.chol, vars, ind=NULL, verbose=0, max.threads=0,seed=NULL,LDL=FALSE)
{

  if(missing(mu))
	  stop('Must specify mean value')

  if(missing(Q) && missing(Q.chol))
	  stop('Must specify a precision matrix or its Cholesky factor')

  if(!missing(ind) && !missing(Q.chol))
	  stop('Cannot provide both cholesky factor and indices.')

  if (!missing(Q.chol) && !is.null(Q.chol)) {
      L = Q.chol
  } else {
      if(LDL){
        L = suppressWarnings(t(as(Cholesky(Q),"Matrix")))
      } else {
        L = chol(Q)
      }
  }

  if(missing(vars)){
    sd <- sqrt(excursions.variances(L))
	}

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

  return(list(a=a,b=b))
}

