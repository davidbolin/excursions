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

excursions <- function(alpha, u, mu, Q, type, n.iter=10000, Q.chol,
                       vars, rho, reo, method='EB', ind, max.size,
                       verbose=0, max.threads=0,seed)
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

  if(missing(mu))
	  stop('Must specify mean value')

  if(missing(Q) && missing(Q.chol))
	  stop('Must specify a precision matrix or its Cholesky factor')

  if(missing(type))
	  stop('Must specify type of excursion set')

  if(qc && missing(rho))
	  stop('rho must be provided if QC is used.')

  if(!missing(ind) && !missing(reo))
	  stop('Either provide a reordering using the reo argument or provied a set of nodes using the ind argument, both cannot be provided')

  if (!missing(Q.chol) && !is.null(Q.chol)) {
      ## make the representation unique (i,j,v)
      Q = Q.chol#private.as.dgTMatrix(Q.chol)
      is.chol = TRUE
  } else {
      ## make the representation unique (i,j,v)
      Q = Q#private.as.dgTMatrix(Q)
      is.chol = FALSE
  }

  if (missing(vars)) {
    L = chol(Q)
    vars <- excursions.variances(L)
  }

  if(verbose)
    cat("Calculate marginals\n")
  marg <- excursions.marginals(type = type, rho = rho,vars = vars,
                               mu = mu, u = u, QC = qc, ind = ind)

  if (missing(max.size)){
	  m.size = length(mu)
  } else {
	  m.size = max.size
  }
  if (!missing(ind)) {
	  indices = rep(FALSE,length(mu))
	  indices[ind] = TRUE
	  if(missing(max.size)){
		  m.size = length(ind)
	  } else {
		  m.size = min(length(ind),m.size)
	  }
  } else {
    indices = rep(TRUE,length(mu))
  }

  if(verbose)
    cat("Calculate permutation\n")
  if(missing(reo)){
    use.camd = !missing(ind) || alpha < 1
    if(qc){
      reo <- excursions.permutation(marg$rho_ng, indices, use.camd = TRUE,alpha,Q)
    } else {
      reo <- excursions.permutation(marg$rho, indices, use.camd = TRUE,alpha,Q)
    }
  }

  if(verbose)
    cat("Calculate limits\n")
  limits <- excursions.setlimits(marg, vars,type,QC=qc,u,mu)

  res <- excursions.call(limits$a,limits$b,reo,Q, is.chol = is.chol,
                         1-alpha, K = n.iter, max.size = m.size,
                         n.threads = max.threads,seed = seed)

  n = length(mu)
  ii = which(res$Pv[1:n] > 0)
  if (length(ii) == 0) i=n+1 else i=min(ii)

  F = Fe  = rep(0,n)
  F[reo] = res$Pv
  Fe[reo] = res$Ev

  ireo = NULL
  ireo[reo] = 1:n

  if(type == '='){
	  F=1-F
	  D = rep(1,n)
	  if(i<n+1) D[reo[i:n]] = 0
  } else {
	  D = rep(0,n)
	  if(i<n+1) D[reo[i:n]] = 1
  }

  G = rep(-1,n)
  if(type == ">"){
    G[mu>u] = 1
  } else if(type == "<") {
    G[mu<u] = 0
  } else if(type == "=" || type == "!") {
    G[mu<u] = 0
    G[mu>u] = 1
  }

  output <- list(F=F, G = G,
                 D=D, rho=marg$rho, Fe=Fe, reo=reo, ireo=ireo, vars=vars,
                 meta=list(calculation="excursions",
                 type=type,
                 level=u,
                 alpha=alpha,
                 n.iter=n.iter,
                 method=method))
  class(output) <- "excurobj"
  output
}
