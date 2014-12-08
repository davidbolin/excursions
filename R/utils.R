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




## calculates variances given either a cholesky factor L in Matrix or SPAM f
## format, or given a precision matrix Q. If Q is provided, the cholesky factor
## is calculated and the variances are then returned in the same ordering as Q
## If L is provided, the variances are returned in the same ordering as L, even
## if L@invpivot exists.
excursions.variances<-function(L,Q)
{

  if(!missing(L) && !is.null(L)){
    ireo = FALSE
    if(is(L,'spam.chol.NgPeyton')){
     L = as(as(as.dgRMatrix.spam(as.spam(L)), "TsparseMatrix"),"dtCMatrix")
   } else {
      if (!is(L, "dtCMatrix"))
       stop("L needs to be in ccs format for now.")
    }
  } else {
    L = chol(private.as.spam(Q))
    ireo = TRUE
    reo = L@invpivot
    L = as(as(as.dgRMatrix.spam(as.spam(L)), "TsparseMatrix"),"dtCMatrix")

  }
  Mp = L@p
  Mi = L@i
  Mv = L@x
  n = dim(L)[1]

  out<- .C("Qinv",Rir = as.integer(Mi), Rjc = as.integer(Mp),
           Rpr = as.double(Mv), variances=double(n), n = as.integer(n))
  if(ireo){
    return(out$variances[reo])
  } else {
    return(out$variances)
  }


}


excursions.marginals <- function(type, rho,vars, mu, u, QC = FALSE)
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


excursions.call <- function(a,b,reo,Q, is.chol = FALSE, lim, K, max.size,n.threads, seed,LDL=FALSE)
{
  if(is.chol == FALSE){
    a.sort = a[reo]
    b.sort = b[reo]
    Q = Q[reo,reo]

    if(LDL){
  		L = suppressWarnings(t(as(Cholesky(Q,perm=FALSE),"Matrix")))
    } else {
      L = suppressWarnings(chol.spam(private.as.spam(Q),pivot = FALSE))
    }

    res = gaussint(Q.chol = L, a = a.sort, b = b.sort, lim = lim,
                                 n.iter = K, max.size = max.size,
                                 max.threads = n.threads, seed = seed)
  } else {
    #assume that everything already is ordered
    res = gaussint(Q = Q, a= a, b = b, lim = lim, n.iter = K,
                    max.size = max.size,
                    max.threads = n.threads,seed = seeed,LDL=LDL)
  }
  return(res)
}

private.as.spam = function(A)
{
  if(is(A,"spam")){
    return(A)
  }
  else if(is(A,"dsyMatrix")){
    return(as.spam.dgCMatrix(as(as.matrix(A),"dgCMatrix")))
  } else {
    return(as.spam.dgCMatrix(as(A,"dgCMatrix")))
  }
}


##
# Distribution function of Gaussian mixture \Sum_k w[k]*N(mu[k],sigma[k]^2)
##
Fmix = function(x,mu,sd,w) sum(w*pnorm(x,mean=mu,sd=sd))

##
# Quantile function of Gaussian mixture
##
Fmix_inv = function(p,mu,sd,w,br=c(-1000,1000))
{
   G = function(x) Fmix(x,mu,sd,w) - p
   return(uniroot(G,br)$root)
}

##
# Function for optimization of interval for mixtures
##
fmix.opt <- function(x,
                     alpha,
                     sd,
                     Q.chol,
                     w,
                     mu,
                     limits,
                     verbose,
                     max.threads,
                     ind)
{
  K = dim(mu)[1]
  n = dim(mu)[2]
  q.a = sapply(seq_len(n),function(i) Fmix_inv(x/2,
                                               mu = mu[,i],
                                               sd = sd[,i],
                                               w=w,
                                               br=limits))
  q.b = sapply(seq_len(n),function(i) Fmix_inv(1-x/2,
                                               mu = mu[,i],
                                               sd = sd[,i],
                                               w=w,
                                               br=limits))



  prob = 0
  stopped = 0

  k.seq = sort(w,decreasing=TRUE,index.return=TRUE)$ix

  ki = 1
  for(k in k.seq){
    ws = 0
    if(ki<K){
      ws = sum(w[k.seq[(ki+1):K]])
    }
    ki = ki+1
    lim = (1-alpha - prob - ws)/w[k]
    p = gaussint(mu = mu[k,],
                 Q.chol=Q.chol[[k]],
                 a=q.a,
                 b=q.b,
                 ind = ind,
                 lim=max(0,lim),
                 max.threads=max.threads)
    if(p$P==0){
      stopped = 1
      break
    } else {
      prob = prob + w[k]*p$P
    }
  }


  if(stopped==1){ #too large alpha
    if(prob == 0){
      val = 10*(1+x)
    } else {
      val = 1+(prob-(1-alpha))^2
    }
  } else { #too small x
    val = (prob-(1-alpha))^2
  }

  if(verbose){
    cat("in optimization: ",x," ", prob, " ", val, "\n")
  }

  return(val)

}


fmix.samp.opt <- function(x, alpha,mu, sd, w, limits, samples)
{
  n = dim(mu)[2]
  q.a = sapply(seq_len(n),function(i) Fmix_inv(x/2,
                                               mu = mu[,i],
                                               sd = sd[,i],
                                               w=w,
                                               br=limits))
  q.b = sapply(seq_len(n),function(i) Fmix_inv(1-x/2,
                                               mu = mu[,i],
                                               sd = sd[,i],
                                               w=w,
                                               br=limits))

  cover = sapply(seq_len(dim(samples)[1]), function(i) (sum(samples[i,]>q.b) + sum(samples[i,]<q.a))==0)

  prob = mean(cover)
  val = (prob-(1-alpha))^2
  cat("in optimization: ",x," ", prob, " ", val, "\n")
  return(val)

}

