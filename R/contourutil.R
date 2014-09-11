## Calculate the contour map function
contourfunction <- function(lp,mu,Q,vars,ind, alpha=1, n.iter=10000,
                                       Q.chol,max.threads=0)
{
	if (!missing(Q.chol) && !is.null(Q.chol)) {
      Q = Q.chol
      is.chol = TRUE
  } else if (!missing(Q) && !is.null(Q)) {
      is.chol = FALSE
  } else {
    stop('Must supply Q or Q.chol')
  }

  if(missing(mu))
    stop('Must supply mu')

  if(missing(lp))
    stop('Must supply level plot')

	lim <- excursions.limits(lp=lp,mu=mu,measure=0,ind=ind)

  if (missing(vars)) {
    if(is.chol) {
      vars <- excursions.variances(Q)
    } else {
      L = chol(Q)
      vars <- excursions.variances(L)
    }
  }
  rho <- contourmap.marginals(mu=mu,vars=vars,lim=lim,ind=ind)

  lim$a <- lim$a - mu
	lim$b <- lim$b - mu

	m.size = length(mu)
  indices = NULL
  if (!missing(ind)) {
	  indices = rep(0,length(mu))
	  indices[ind] = 1
		m.size = length(ind)
  }

  use.camd = !missing(ind) || alpha < 1
  reo <- excursions.permutation(rho = rho, ind = indices,
                                use.camd = use.camd,alpha = alpha,Q = Q)

  res <- excursions.call(lim$a,lim$b,reo,Q, is.chol = is.chol,
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

	D = rep(0,n)
	if(i<n+1) D[reo[i:n]] = 1

  return(list(F=F, Fe=Fe, D=D, rho=rho))
}

## Calculate marginal probabilities P(lim$a < X < lim$b) for
## when X is N(mu,vars)-distributed
contourmap.marginals <- function(mu,vars,lim,ind)
{
  if(!missing(ind) && !is.null(ind)){
	  marg = rep(0,length(mu))
	  marg[ind] = pnorm(lim$b[ind], mu[ind], sqrt(vars[ind])) -
	              pnorm(lim$a[ind], mu[ind], sqrt(vars[ind]))
	} else {
	  marg = pnorm(lim$b, mu, sqrt(vars)) - pnorm(lim$a, mu, sqrt(vars))
	}
	return(marg)
}

## Create a levelplot with given levels/number of levels
excursions.levelplot <- function(mu,n.levels,ind,levels,equal.area=FALSE)
{
	n = length(mu)
	if(missing(ind)) ind = 1:n
	x.mean = rep(-Inf,n)
	x.mean[ind] = mu[ind]
	r = range(mu[ind])
	if(missing(levels)){
		if(missing(n.levels)){
			stop('Must specify levels or number of levels')
		}
		if(equal.area){
			levels <- as.vector(quantile(mu[ind],
								(1:n.levels)/(n.levels+1),type=4))
		} else {
			levels <- seq(from=r[1],to=r[2],
						  length.out = (n.levels+2))[2:(n.levels+1)]
		}
	} else {
		if(!missing(n.levels)){
			if(n.levels != length(levels)){
			  warning('n.levels is not equal to the length of levels')
				n.levels = length(levels)
			}
		} else {
			n.levels = length(levels)
		}
	}

	u.e = NULL
	E = vector("list",n.levels+1)
	l1 = c(r[1],levels,r[2])
	for(i in 1:(n.levels+1)) u.e[i] = (l1[i]+l1[i+1])/2

	for(i in 1:(n.levels)){
		E[[i]] = which((l1[i] <= x.mean) & (x.mean < l1[i+1]))
	}
	E[[n.levels+1]] = which((l1[n.levels+1] <= x.mean))
	map = rep(0,n)
	for(i in 1:(n.levels+1)) map[E[[i]]] = u.e[i]

	return(list(u = levels, n.levels = n.levels, u.e = u.e, E=E,map=map))
}

## Create a P-optimal levelplot.
## The function will take A LOT of time to run if use.marginals=FALSE.
excursions.opt.levelplot <- function(mu,vars,Q,n.levels, measure=2, use.marginals=TRUE, ind)
{
	if( (measure != 1) && (measure != 2) && (measure != 0))
		stop('only measure 0, 1, or 2 allowed')

	if(missing(ind)) ind=1:length(mu)

	r = range(mu[ind])

	start.v <- seq(from=r[1],to=r[2],length.out = (n.levels+2))[2:(n.levels+1)]
  Q.chol = chol(Q)
	u = optim(start.v,excursions.lim.func,mu=mu,vars=vars,Q.chol=Q.chol,
			  measure=measure,use.marginals = use.marginals,ind=ind,Q=Q)

	lp = excursions.levelplot(mu,levels = u$par,ind=ind)
	if(measure==2){
		if(use.marginals){
			lp$P2bound = -u$value
		} else {
			lp$P2 = -u$value
		}
	} else if(measure==1){
		if(use.marginals){
			lp$P1bound = -u$value
		} else {
			lp$P1 = -u$value
		}
	} else if(measure==0){
		if(use.marginals){
			lp$P0bound = -u$value
		} else {
			lp$P0 = -u$value
		}
	}
	return(lp)
}

## Internal function for optimization of P-optimal contour map
excursions.lim.func <- function(u, mu, vars, Q.chol, Q, measure,
                                use.marginals, ind=ind)
{
	lp = excursions.levelplot(mu,levels = u,ind=ind)
	if( min(u)<range(mu[ind])[1] | max(u) > range(mu[ind])[2]){
	  # levels should be in (min(mu),max(mu))
	  val = 0
	} else if (max(sort(u,index.return=TRUE)$ix - seq_len(length(u)))>0) {
	  # levels should be sorted
	  val = 0
	} else {
	  v = TRUE;
	  for(i in 1:(length(u)-1)){
	    v = v & sum(mu[ind]<u[i+1] & mu[ind] > u[i])>0
	  }
	  if(v == FALSE) {
	    # all sets E should be non-empty
	    val = 0
	  } else {
  		if(use.marginals){
	  		limits = excursions.limits(lp,mu,measure=measure)
		  	val = -min(contourmap.marginals(mu,vars,limits,ind)[ind])
		  } else {
			  val = -Pmeasure(lp,mu=mu,Q=Q,Q.chol=Q.chol,type=measure,
			                ind=ind,vars=vars)
	      cat(u, ': ', -val, '\n')
		  }
		}
	}
	return(val)
}

## Function that calculates the P measure for a given contour map.
Pmeasure <- function(lp,mu,Q,Q.chol, ind=NULL,type,vars=vars)
{
  if(type==0){
    res <- contourfunction(lp=lp,mu=mu,Q=Q,vars=vars,ind=ind)
    p = mean(res$F)
  } else {
    limits = excursions.limits(lp=lp,mu=mu,measure=type)
	  res = gaussint(mu = mu, Q=Q, Q.chol = Q.chol, a=limits$a,
	                b=limits$b,ind=ind,reo=TRUE)
	  p = res$P[1]
  }
	return(p)
}

## Set integration limits for a given measure.
excursions.limits <- function(lp,mu,measure,ind)
{
	n = length(mu)
	n.l = length(lp$u)
	if(measure == 1 & n.l<2){
	  stop('P1 measure only makes sense if number of contours is larger than 1')
	}

	a = rep(-Inf,n)
	b = rep(Inf,n)

	if(measure==2){
		b[lp$E[[1]]] = lp$u.e[2]
		a[lp$E[[n.l+1]]] = lp$u.e[n.l]

		if(n.l>1){
			for(i in 2:n.l){
				a[lp$E[[i]]] = lp$u.e[i-1]
				b[lp$E[[i]]] = lp$u.e[i+1]
			}
		}
	} else if(measure ==1){
		if(n.l>1){
			b[lp$E[[1]]] = lp$u[2]
			a[lp$E[[n.l+1]]] = lp$u[n.l-1]
		}
		if(n.l>2){
			b[lp$E[[2]]] = lp$u[3]
			a[lp$E[[n.l]]] = lp$u[n.l-2]
		}
		if(n.l>3){
			for(i in 3:(n.l-1)){
				a[lp$E[[i]]] = lp$u[i-2]
				b[lp$E[[i]]] = lp$u[i+1]
			}
		}
	} else if(measure==0){

		b[lp$E[[1]]] = lp$u[1]
		a[lp$E[[n.l+1]]] = lp$u[n.l]

		if(n.l>1){
			for(i in 2:n.l){
				a[lp$E[[i]]] = lp$u[i-1]
				b[lp$E[[i]]] = lp$u[i]
			}
		}
	} else {
		stop('Measure must be 0, 1, or 2')
	}
	return(list(a=a,b=b))
}
