
excursions.integration <- function(L, a, b, lim = 0, n.iter = 10000, max.size, n.threads=0,seed)
{

  if (!is(L, "dtCMatrix"))
      stop("L needs to be in ccs format for now.")
    
  n = dim(L)[1]
  Mp = L@p
  Mi = L@i
  Mv = L@x
  
  if(missing(max.size))
    max.size = n
  

  if(missing(seed)){
    seed_provided = 0
    seed.in = as.integer(rep(0,6))
  } else {
    seed_provided = 1
    seed.in = seed
  }

  Pv = rep(0,n)
  Ev = rep(0,n)
  
  opts = c(n,n.iter,max.size,n.threads,seed_provided)

  out <- .C("shapeInt", Mp = as.integer(Mp), Mi = as.integer(Mi), 
              Mv = as.double(Mv), a = as.double(a), b = as.double(b), 
              opts = as.integer(opts), lim = as.double(lim),
              Pv = as.double(Pv), Ev = as.double(Ev),seed_in=seed.in)

  return(list(P = out$Pv, E = out$Ev)) 
}

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
  if(!missing(ind)){
    rl$rho[ind==0] = -1
    if(QC)
     rl$rho_ng[ind==0] = -1
  }
  return(rl)
}

	
excursions.permutation <- function(rho, ind, use.camd = TRUE,alpha,Q)
{
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
    a.sort[a.sort==Inf]  = .Machine$double.xmax
    b.sort[b.sort==Inf]  = .Machine$double.xmax
    a.sort[a.sort==-Inf] = -.Machine$double.xmax
    b.sort[b.sort==-Inf] = -.Machine$double.xmax

    res = excursions.integration(L, a.sort, b.sort, lim, K, max.size, n.threads, seed)
  } else {
    #assume that everything already is ordered
    a[a==Inf]  = .Machine$double.xmax
    b[b==Inf]  = .Machine$double.xmax
    a[a==-Inf] = -.Machine$double.xmax
    b[b==-Inf] = -.Machine$double.xmax
    res = excursions.integration(Q, a, b, lim, K, max.size, n.threads,seed)
  }
  return(res)
}
	
	