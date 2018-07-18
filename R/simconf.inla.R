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

#' Simultaneous confidence regions for latent Gaussian models
#'
#' \code{simconf.inla} is used for calculating simultaneous confidence regions
#' for latent Gaussian models estimated using \code{INLA}.
#'
#' @param result.inla Result object from INLA call.
#' @param stack The stack object used in the INLA call.
#' @param name The name of the component for which to do the calculation. This
#' argument should only be used if a stack object is not provided, use the tag
#' argument otherwise.
#' @param tag The tag of the component in the stack for which to do the calculation.
#' This argument should only be used if a stack object is provided, use the name
#' argument otherwise.
#' @param ind If only a part of a component should be used in the calculations, this
#' argument specifies the indices for that part.
#' @param alpha Error probability for the region.
#' @param method Method for handeling the latent Gaussian structure:
#' \itemize{
#' \item{'EB' }{Empirical Bayes (Gaussian approximation of posterior).}
#' \item{'NI' }{Numerical integration (Calculation based on the Gaussian mixture
#' approximation of the posterior, as calculated by INLA).}}
#' @param n.iter Number or iterations in the MC sampler that is used for approximating
#' probabilities. The default value is 10000.
#' @param verbose Set to TRUE for verbose mode (optional).
#' @param link Transform output to the scale of the data using the link function as defined in
#' the model estimated with INLA (default FALSE).
#' @param max.threads Decides the number of threads the program can use. Set to 0 for
#' using the maximum number of threads allowed by the system (default).
#' @param seed Random seed (optional).
#' @param inla.sample Set to TRUE if inla.posterior.sample should be used for the MC
#' integration.
#'
#' @return An object of class "excurobj" with elements
#' \item{a }{The lower bound.}
#' \item{b }{The upper bound.}
#' \item{a.marginal }{The lower bound for pointwise confidence bands.}
#' \item{b.marginal }{The upper bound for pointwise confidence bands.}
#' @export
#' @details See \code{\link{simconf}} for details.
#' @note This function requires the \code{INLA} package, which is not a CRAN package.
#' See \url{http://www.r-inla.org/download} for easy installation instructions.
#' @author David Bolin \email{davidbolin@@gmail.com}
#' @references Bolin et al. (2015) \emph{Statistical prediction of global sea level
#' from global temperature}, Statistica Sinica, Vol 25, pp 351-367.
#' @seealso \code{\link{simconf}}, \code{\link{simconf.mc}}, \code{\link{simconf.mixture}}
#' @examples
#' \donttest{
#' if (require.nowarnings("INLA")) {
#' n <- 10
#' x <- seq(0, 6, length.out=n)
#' y <- sin(x) + rnorm(n)
#' mu <- 1:n
#' result <- inla(y ~ 1 + f(mu, model='rw2'),
#'                data=list(y=y, mu=mu), verbose=FALSE,
#'                control.compute = list(config=TRUE),
#'                num.threads = 1)
#'
#' res <- simconf.inla(result, name='mu', alpha = 0.05, max.threads = 1)
#'
#' plot(result$summary.random$mu$mean,ylim=c(-2,2))
#' lines(res$a)
#' lines(res$b)
#' lines(res$a.marginal,col="2")
#' lines(res$b.marginal,col="2")
#' }}

simconf.inla <- function(result.inla,
                         stack,
                         name=NULL,
                         tag=NULL,
                         ind=NULL,
                         alpha,
                         method="NI",
                         n.iter=10000,
                         verbose=0,
                         link = FALSE,
                         max.threads=0,
                         seed=NULL,
                         inla.sample=TRUE)
{
  if (!requireNamespace("INLA", quietly=TRUE))
    stop('This function requires the INLA package (see www.r-inla.org/download)')
  if(missing(result.inla))
    stop('Must supply INLA result object')

  if(missing(alpha))
    stop('Must supply error probability alpha')

  if(result.inla$.args$control.compute$config==FALSE)
    stop('INLA result must be calculated using control.compute$config=TRUE')

  n = length(result.inla$misc$configs$config[[1]]$mean)

  if(!missing(ind))
    ind <- private.as.vector(ind)


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

  links = NULL
  if(link){
    links = result.inla$misc$linkfunctions$names[
                                            result.inla$misc$linkfunctions$link]
    links = links[ind]
  }

  n.theta = result.inla$misc$configs$nconfig

  if(method=="EB") {
    for(i in 1:n.theta){
      config = private.get.config(result.inla,i)
      if(config$lp == 0)
        break
    }
    res <- simconf(alpha = alpha, mu = config$mu, Q = config$Q,
                   vars = config$vars, n.iter=n.iter, ind=ind,
                   verbose=verbose, max.threads=max.threads, seed=seed)
    res$meta$call = match.call()
    return(private.simconf.link(res,links,link))


  } else if(method=="NI") {
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
      s = suppressWarnings(INLA::inla.posterior.sample(n.iter,result.inla))
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
      limits = c(-1000,1000)
      a.marg = sapply(seq_len(length(ind)), function(i) Fmix_inv(p = alpha/2,
                                                   mu = mu.m[,i], sd = sd.m[,i],
                                                   w = w, br = limits))

      b.marg = sapply(seq_len(length(ind)), function(i) Fmix_inv(p = 1-alpha/2,
                                                   mu = mu.m[,i], sd = sd.m[,i],
                                                   w = w, br = limits))
      while(min(a.marg) == limits[1] || max(b.marg) == limits[2])
      {
        limits = 2*limits
        a.marg = sapply(seq_len(length(ind)), function(i) Fmix_inv(p = alpha/2,
                                                   mu = mu.m[,i], sd = sd.m[,i],
                                                   w = w, br = limits))

        b.marg = sapply(seq_len(length(ind)), function(i) Fmix_inv(p= -alpha/2,
                                                   mu = mu.m[,i], sd = sd.m[,i],
                                                   w = w, br = limits))
      }

      r.o = optimize(fmix.samp.opt,interval = c(0,alpha), mu=mu.m, alpha=alpha,
                     sd=sd.m, w=w, limits = limits,samples=samp)

      a = sapply(seq_len(length(ind)), function(i) Fmix_inv(p = r.o$minimum/2,
                                              mu = mu.m[,i], sd = sd.m[,i],
                                              w = w, br = limits))

      b = sapply(seq_len(length(ind)), function(i) Fmix_inv(p = 1-r.o$minimum/2,
                                              mu = mu.m[,i], sd = sd.m[,i],
                                              w = w, br = limits))

      res = list(a = a, b = b, a.marginal = a.marg, b.marginal = b.marg,
                 mean = mu, vars = vars)

      res$meta = list(calculation="simconf",
                      alpha=alpha,
                      n.iter=n.iter,
                      ind=ind,
                      call = match.call())
      class(res) <- "excurobj"

      return(private.simconf.link(res,links,link))

    } else {
      res = simconf.mixture(alpha = alpha, mu = mu, Q = Q, vars = vars,
                            w = w, n.iter=n.iter, ind=ind, verbose=verbose,
                            max.threads=max.threads, seed=seed)
      res$meta$call = match.call()
    return(private.simconf.link(res,links,link))
    }
  } else {
    stop("Method must be EB or NI")
  }
}
