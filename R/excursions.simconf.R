## excursions.simconf.R
##
##   Copyright (C) 2012, 2013 David Bolin, Finn Lindgren
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

excursions.simconf <- function(alpha,mu,Q,n.iter=10000,Q.chol,vars,ind,verbose=0,max.threads=0){

# Extract input data and setup parameters
if(missing(mu))
	stop('Must specify mean value')

if(missing(Q) && missing(Q.chol))
	stop('Must specify a precision matrix or its Cholesky factor')

if(!missing(ind) && !missing(Q.chol))
	stop('Cannot provide both cholesky factor and indices.')

if(missing(vars)){
	if(require("INLA")){
		sd = sqrt(diag(INLA::inla.qinv(Q)))
	} else {
		sd = sqrt(diag(chol2inv(chol(Q))))
	}
} else {
	sd = sqrt(vars)
}

#setup function for optmization
f.opt <- function(x,alpha,sd,Q,Q.chol,ind){
	q = qnorm(x)*sd;
	prob = excursions.int(mu = rep(0,Q@Dim[1]), Q=Q, a=-q, b=q,
	 					  Q.chol=Q.chol, ind=ind, alpha=alpha,
	 					  max.threads=max.threads)

	if(prob$P == 0){
		return(1)
	} else {
		return(prob$P)
	}

}

r.o = optimize(f.opt,interval = c(0,1),alpha=alpha,sd=sd,Q=Q)

a = mu-qnorm(r.o$minimum)*sd
b = mu+qnorm(r.o$minimum)*sd

return(list(a=a,b=b))

}

