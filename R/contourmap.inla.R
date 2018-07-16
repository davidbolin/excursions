## excursions.inla.R
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


contourmap.inla <- function(result.inla,
                            stack,
                            name=NULL,
                            tag=NULL,
                            method="QC",
                            n.levels,
                            type = c("standard",
                                     "pretty",
                                     "equalarea"),
                            compute = list(F=TRUE,measures=NULL),
                            alpha,
                            F.limit,
                            n.iter=10000,
                            verbose=FALSE,
                            max.threads=0,
                            seed=NULL,
                            ind,...)
{
  if (!requireNamespace("INLA", quietly=TRUE))
    stop('This function requires the INLA package (see www.r-inla.org/download)')
  if(missing(result.inla))
    stop('Must supply INLA result object')

  type <- match.arg(type)

  if(missing(alpha) || is.null(alpha)){
    alpha = 0.1
  }
  if(missing(F.limit)) {
    F.limit = 0.99
  } else {
    F.limit = max(alpha,F.limit)
  }
  if(missing(n.levels) || is.null(n.levels)){
    stop("Must supply n.levels")
  }
  if(!is.null(compute$measures))
    measure <- match.arg(compute$measures,
                         c("P0", "P1", "P2","P0-bound","P1-bound","P2-bound"),
                         several.ok=TRUE)
  if(compute$F){ #compute P0 measure if F is computed anyway
    if(("P0" %in% measure)==FALSE){
      measure = c(measure,"P0")
    }
  }
  if(method != 'EB' && method != 'QC' )
    stop("Currently only EB and QC methods are implemented")

  #compute indices, here ind will contain the indices that are used to extract the
  #relevant part from the configs, ind.int is the index vector for extracting marginal
  #distributions for random effects, and indices is a logical version of ind
  ind.stack <- inla.output.indices(result.inla, name=name, stack=stack, tag=tag)
  n.out <- length(ind.stack)
  ind.int <- seq_len(n.out)
  if(!missing(ind) && !is.null(ind)){
    ind.int <- ind.int[ind]
    ind.stack <- ind.stack[ind]
  }
  ind = ind.stack

  for(i in 1:result.inla$misc$configs$nconfig){
    config = private.get.config(result.inla,i)
    if(config$lp == 0)
      break
  }
  indices = rep(FALSE,length(config$mu))
  indices[ind] = TRUE

  cm <- contourmap(mu=config$mu,Q = config$Q,
                   ind=ind,
                   compute = list(F=FALSE,measures=NULL),
                   n.levels = n.levels,...)

  #Are we interested in a random effect?
  random.effect = FALSE
  if(!missing(name)&&!is.null(name)&&name!="APredictor"&&name!="Predictor"){
    random.effect = TRUE
    if(is.null(result.inla$marginals.random))
      stop('INLA result must be calculated using return.marginals.random=TRUE if P measures are to be calculated for a random effect of the model')
  }

  if(!random.effect && is.null(result.inla$marginals.linear.predictor))
    stop('INLA result must be calculated using return.marginals.linear.predictor=TRUE if P measures are to be calculated for the linear predictor.')


  #compute measures
  F.computed = FALSE
  if(!is.null(measure)){
    rho <- matrix(0,length(config$mu),2)
    for( i in 1:length(measure)){
      #compute limits
      limits = excursions.limits(lp=cm,mu=config$mu,measure=measure[i])
      #if QC, update limits
      if(method=="QC"){
        if(random.effect) {
          rho.ind <- sapply(1:length(ind), function(j) inla.get.marginal.int(ind.int[j],
                                                                             a = limits$a[ind[j]],
                                                                             b = limits$b[ind[j]],
                                                                             result = result.inla,
                                                                             effect.name=name))
        } else {
          rho.ind <- sapply(1:length(ind), function(j) inla.get.marginal.int(ind[j],
                                                                             a = limits$a[ind[j]],
                                                                             b = limits$b[ind[j]],
                                                                             result = result.inla))
        }
        rho[ind,] = t(rho.ind)
        limits$a <- config$mu + sqrt(config$vars)*qnorm(pmin(pmax(rho[,1],0),1))
        limits$b <- config$mu + sqrt(config$vars)*qnorm(pmin(pmax(rho[,2],0),1))
      } else {
        rho = NULL
      }

      if(measure[i]=="P1") {
        if(n.levels>1){
          if(verbose) cat('Calculating P1-measure\n')
          tmp = gaussint(mu = config$mu, Q=config$Q, a=limits$a,
                         b=limits$b,ind=indices,use.reordering="limits",
                         n.iter=n.iter,seed=seed)
          cm$P1 <- tmp$P[1]
          cm$P1.error <- tmp$E[1]
        } else {
          cm$P1 = 1
          cm$P1.error <- 0
        }
      } else if(measure[i] == "P2") {
        if(verbose) cat('Calculating P2-measure\n')
        tmp = gaussint(mu = config$mu, Q=config$Q, a=limits$a,
                       b=limits$b,ind=indices,use.reordering="limits",
                       n.iter=n.iter,seed=seed)
        cm$P2 <- tmp$P
        cm$P2.error <- tmp$E
      } else if (measure[i] == "P0") {
        if(verbose) cat('Calculating P0-measure and contour map function\n')
        p <- contourfunction(lp=cm, mu=config$mu,Q=config$Q ,vars=config$vars,
                             ind = ind,alpha=alpha, F.limit = F.limit,
                             rho = rho,n.iter=n.iter,max.threads=max.threads,
                             seed=seed,verbose=verbose)
        cm$P0 = mean(p$F[ind])
        cm$F = p$F
        cm$E = p$E
        cm$M = p$M
        cm$rho = p$rho
        cm$meta$F.limit=F.limit
        cm$meta$alpha=alpha
        F.computed = TRUE
      } else{ #bounds
        if(method!="QC"){
          rho = matrix(0,length(config$mu),2)
          rho[ind,1] = pnorm(limits$a[ind], config$mu[ind], sqrt(config$vars[ind]))
          rho[ind,2] = pnorm(limits$b[ind], config$mu[ind], sqrt(config$vars[ind]))
        }
        if(measure[i] == "P0-bound"){
          cm$P0.bound <- mean(rho[ind,2]-rho[ind,1])
        } else if(measure[i] == "P1-bound"){
          cm$P1.bound <- min(rho[ind,2]-rho[ind,1])
        } else if(measure[i] == "P2-bound"){
          cm$P2.bound <- min(rho[ind,2]-rho[ind,1])
        }
      }
    }
  }
  #Fix output
  set.out = rep(NA,n.out)
  if(F.computed){
    set.out[ind.int] = cm$E[ind];
    cm$E = set.out
  }

  if(!is.null(cm$G)) {
    set.out[ind.int] = cm$G[ind];
    cm$G = set.out
  }

  if(!is.null(cm$F)){
    F0 = cm$F[ind]; F0[is.na(F0)] = 0
    cm$P0 = mean(F0)
    set.out[ind.int] = cm$F[ind];
    cm$F = set.out
  }
  cm$meta$F.computed = F.computed

  if(!is.null(cm$M))
    set.out[ind.int] = cm$M[ind]; cm$M = set.out

  if(!is.null(cm$rho))
    set.out[ind.int] = cm$rho[ind]; cm$rho = set.out

  if(!is.null(cm$map))
    set.out[ind.int] = cm$map[ind]; cm$map = set.out

  cm$meta$ind <- ind.int
  cm$meta$call = match.call()
  return(cm)
}
