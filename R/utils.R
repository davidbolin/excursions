## utils.R 
##
##   Copyright (C) 2013 David Bolin, Finn Lindgren
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


excursions.variances<-function(L)
{
  if (!is(L, "dtCMatrix"))
      stop("L needs to be in ccs format for now.")
  Mp = L@p
  Mi = L@i
  Mv = L@x 
  n = dim(L)[1]

  out<- .C("Qinv",Rir = as.integer(L@i), Rjc = as.integer(L@p), 
           Rpr = as.double(L@x), variances=double(n), n = as.integer(n))

  return(out$variances)
}

excursions.marginals <- function(type, rho,vars, mu, u, QC = FALSE, ind)	
{  
  rl = list()
  if(type == "=" || type == "!="){
    if(QC){
      rl$rho_ngu = rho
      rl$rho_l = pnorm(mu-u,sd=sqrt(vars), lower.tail=FALSE)
			rl$rho_u = 1-rl$rho_l
			rl$rho = pmax(rl$rho_u,rl$rho_l)
			rl$rho_ng = pmax(rl$rho_ngu,1-rl$rho_ngu)
    } else {
      if(!missing(rho)){
        rl$rho_u = rho  
      } else {
        rl$rho_u = 1 - pnorm(mu-u,sd=sqrt(vars), lower.tail=FALSE)	
		  }
		  rl$rho_l = 1-rl$rho_u
		  rl$rho = pmax(rl$rho_u,rl$rho_l)
    }
  } else {
    if(QC) {
      rl$rho_ng = rho
      if(type == ">"){
			  rl$rho = pnorm(mu-u,sd=sqrt(vars))
			} else {
				rl$rho = pnorm(mu-u,sd=sqrt(vars), lower.tail=FALSE)
			}
    } else {
      if(missing(rho)){
       if(type == ">"){
  				rl$rho = pnorm(mu-u,sd=sqrt(vars))
       } else {
  				rl$rho = pnorm(mu-u,sd=sqrt(vars), lower.tail=FALSE)
        }
      } 
    }
  } 
  return(rl)
}

	
excursions.permutation <- function(rho, ind, use.camd = TRUE,alpha,Q)
{
  if(!missing(ind) && !is.null(ind))
    rho[!ind] = -1
  
  n = length(rho)
  v.s = sort(rho,index.return=TRUE)  
  reo = v.s$ix
  rho_sort = v.s$x
  ireo = integer(n)
  ireo[reo] = 1:n
  if(use.camd){
    k=0
    i = n-1
    #add nodes to lower bound
    cindr = cind = rep(0,n)
    while(rho_sort[i]> 1 - alpha && i>0){
      cindr[i] = k 
			i = i-1
			k=k+1
    }
    if(i>0){
      #reorder nodes below the lower bound for sparsity
      while (i>0) {
				cindr[i] = k;
				i = i-1
			}
      #change back to original ordering
      for (i in 1:n) {
				cind[i] = k - cindr[ireo[i]]
			}
			#call CAMD
			out <- .C("reordering",nin = as.integer(n), Mp = as.integer(Q@p), 
		                        Mi = as.integer(Q@i), reo = as.integer(reo), 
		                        cind = as.integer(cind))
		  reo = out$reo+1  
    } 
  }
  return(reo)
}
	
	
excursions.setlimits <- function(marg, vars,type,QC,u,mu)
{
	if(QC){
		if (type=="<") { 
		  uv = sqrt(vars)*qnorm(marg$rho_ng)	  
		} else if (type==">"){
			uv = sqrt(vars)*qnorm(marg$rho_ng,lower.tail=FALSE)
		} else if (type=="=" || type == "!="){
			uv = sqrt(vars)*qnorm(marg$rho_ngu,lower.tail=FALSE)
		}
	} else {
		uv = u-mu
	}	
	if (type == "=" || type == "!=") {
	  if(QC){
	    a = b = uv
	    a[marg$rho_ngu<=0.5] = -Inf
	    b[marg$rho_ngu>0.5] = Inf
		} else {
		  a = b = uv
		  a[marg$rho_u<=0.5] = -Inf
      b[marg$rho_u>0.5] = Inf
		}
	} else if (type == ">") {
	  a = uv
	  b = rep(Inf,length(mu))
	} else if (type == "<"){
    a = rep(-Inf,length(mu))
    b = uv
  } 
  return(list(a=a,b=b))
}

	
excursions.call <- function(a,b,reo,Q, is.chol = FALSE, lim, K, max.size,n.threads, seed = seed)
{
  if(is.chol == FALSE){
    a.sort = a[reo]
    b.sort = b[reo]
  
    #calculate cholesky here
    L = chol(Q[reo,reo])
    
    res = gaussint(Q.chol = L, a = a.sort, b = b.sort, lim = lim,
                                 n.iter = K, max.size = max.size, 
                                 max.threads = n.threads, seed = seed)
  } else {
    #assume that everything already is ordered
    res = gaussint(Q = Q, a= a, b = b, lim = lim, n.iter = K,
                                 max.size = max.size, 
                                 max.threads = n.threads,seed = seeed)
  }
  return(res)
}
	
