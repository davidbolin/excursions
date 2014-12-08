## simconf.inla.R
##
##   Copyright (C) 2014 David Bolin, Finn Lindgren
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

simconf.inla <- function(result.inla,
                         stack,
                         name=NULL,
                         tag=NULL,
                         ind=NULL,
                         alpha=1,
                         u,
                         n.iter=10000,
                         verbose=0,
                         limits = c(-1000,1000),
                         max.threads=0,
                         seed=NULL,
                         inla.sample=TRUE)
{
  if (!require("INLA"))
    stop('This function requires the INLA package (see www.r-inla.org/download)')
  if(missing(result.inla))
    stop('Must supply INLA result object')

  if(result.inla$.args$control.compute$config==FALSE)
    stop('INLA result must be calculated using control.compute$config=TRUE')

  n = length(result.inla$misc$configs$config[[1]]$mean)

  #Get indices for the component of interest in the configs
  ind.stack <- inla.output.indices(result.inla, name=name, stack=stack, tag=tag)
  n.out = length(ind.stack)
  #Index vector for the nodes in the component of interest
  ind.int <- seq_len(n.out)

  #ind is assumed to contain indices within the component of interest
  if(!missing(ind) && !is.null(ind)){
    ind.int <- ind.int[ind]
    ind.stack <- ind.stack[ind]
  }
  ind = ind.stack

  n.theta = result.inla$misc$configs$nconfig
  mu <- Q <- vars <- list()
  w <- rep(0,n.theta)
  for(i in 1:n.theta){
    config = private.get.config(result.inla,i)
    mu[[i]] = config$mu
    Q[[i]] = config$Q
    vars[[i]] = config$vars
    w[i] = config$lp
  }
  w = exp(w)/sum(exp(w))


  if(inla.sample){
    s = suppressWarnings(inla.posterior.sample(n.iter,result.inla))
    samp <- matrix(0,n.iter,length(ind))

    for(i in seq_len(n.iter)){
      samp[i,] <- s[[i]]$latent[ind]
    }

    mu.m <- matrix(0,n.theta,length(ind))
    sd.m <- matrix(0,n.theta,length(ind))
    for(k in seq_len(n.theta)){
      mu.m[k,] = mu[[k]][ind]
      sd.m[k,] = sqrt(vars[[k]][ind])
    }

    a.marg = sapply(seq_len(length(ind)), function(i) Fmix_inv(p = alpha/2,
                                                   mu = mu.m[,i],
                                                   sd = sd.m[,i],
                                                   w = w,
                                                   br = limits))

    b.marg = sapply(seq_len(length(ind)), function(i) Fmix_inv(p = 1-alpha/2,
                                                   mu = mu.m[,i],
                                                   sd = sd.m[,i],
                                                   w = w,
                                                   br = limits))

    r.o = optimize(fmix.samp.opt,interval = c(0,alpha),
                   mu=mu.m, alpha=alpha,
                   sd=sd.m, w=w, limits = limits,samples=samp)


    a = sapply(seq_len(length(ind)), function(i) Fmix_inv(p = r.o$minimum/2,
                                              mu = mu.m[,i],
                                              sd = sd.m[,i],
                                              w = w,
                                              br = limits))

    b = sapply(seq_len(length(ind)), function(i) Fmix_inv(p = 1-r.o$minimum/2,
                                                          mu = mu.m[,i],
                                                          sd = sd.m[,i],
                                                          w = w,
                                                          br = limits))

  return(list(a = a,
              b = b,
              a.marginal = a.marg,
              b.marginal = b.marg))

  } else {
    return(simconf.mixture(alpha = alpha,
                           mu = mu,
                           Q = Q,
                           vars = vars,
                           w = w,
                           n.iter=n.iter,
                           limits = limits,
                           ind=ind,
                           verbose=verbose,
                           max.threads=max.threads,
                           seed=seed))
  }
}

