## excursions.inla.R
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


#' Excursion sets and contour credible regions for latent Gaussian models
#'
#' Excursion sets and contour credible regions for latent Gaussian models calculated
#' using the INLA method.
#'
#' @param result.inla Result object from INLA call.
#' @param stack The stack object used in the INLA call.
#' @param name The name of the component for which to do the calculation. This argument
#' should only be used if a stack object is not provided, use the tag argument otherwise.
#' @param tag The tag of the component in the stack for which to do the calculation.
#' This argument should only be used if a stack object is provided, use the name argument
#' otherwise.
#' @param ind If only a part of a component should be used in the calculations, this
#' argument specifies the indices for that part.
#' @param method Method for handeling the latent Gaussian structure:
#' \itemize{
#' \item{'EB' }{Empirical Bayes}
#' \item{'QC' }{Quantile correction}
#' \item{'NI' }{Numerical integration}
#' \item{'NIQC' }{Numerical integration with quantile correction}
#' \item{'iNIQC' }{Improved integration with quantle correction}}
#' @param alpha Error probability for the excursion set of interest. The default
#' value is 1.
#' @param F.limit Error probability for when to stop the calculation of the excursion
#' function. The default value is \code{alpha}, and the value cannot be smaller than
#' \code{alpha}. A smaller value of \code{F.limit} results in asmaller compontation time.
#' @param u Excursion or contour level.
#' @param u.link If u.link is TRUE, \code{u} is assumed to be in the scale of the
#' data and is then transformed to the scale of the linear predictor (default FALSE).
#' @param type Type of region:
#'  \itemize{
#'      \item{'>' }{positive excursions}
#'      \item{'<' }{negative excursions}
#'      \item{'!=' }{contour avoiding function}
#'      \item{'=' }{contour credibility function}}
#' @param n.iter Number or iterations in the MC sampler that is used for approximating
#' probabilities. The default value is 10000.
#' @param verbose Set to TRUE for verbose mode (optional).
#' @param max.threads Decides the number of threads the program can use. Set to 0 for using
#' the maximum number of threads allowed by the system (default).
#' @param seed Random seed (optional).
#'
#' @return \code{excursions.inla} returns an object of class "excurobj". This is a list that contains the following arguments:
#' \item{E }{Excursion set, contour credible region, or contour avoiding set}
#' \item{F }{The excursion function corresponding to the set \code{E} calculated
#' for values up to \code{F.limit}}
#' \item{G }{ Contour map set. \eqn{G=1} for all nodes where the \eqn{mu > u}.}
#' \item{M }{ Contour avoiding set. \eqn{M=-1} for all non-significant nodes. \eqn{M=0} for nodes where the process is significantly below \code{u} and \eqn{M=1} for all nodes where the field is significantly above \code{u}. Which values that should be present depends on what type of set that is calculated.}
#' \item{rho }{Marginal excursion probabilities}
#' \item{mean }{Posterior mean}
#' \item{vars }{Marginal variances}
#' \item{meta }{A list containing various information about the calculation.}
#' @export
#' @details The different methods for handling the latent Gaussian structure are listed 
#' in order of accuracy and computational cost. The \code{EB} method is the simplest and is based on a Gaussian 
#' approximation of the posterior of the quantity of interest. The \code{QC} method uses the
#' same Gaussian approximation but improves the accuracy by modifying the limits in the 
#' integrals that are computed in order to find the region. The other three methods are
#' intended for Bayesian models where the posterior distribution for the quantity of 
#' interest is obtained by integrating over the parameters in the model. The \code{NI} 
#' method approximates this integration in the same way as is done in INLA, and the
#' \code{NIQC} and \code{iNIQC} methods combine this apprximation with the QC method
#' for improved accuracy. 
#' 
#' If the main purpose of the analysis is to construct excursion or contour sets for low
#' values of \code{alpha}, we recommend using \code{QC} for problems with Gaussian
#' likelihoods and \code{NIQC} for problems with non-Gaussian likelihoods. The reason for 
#' this is that the more accurate methods also have higher computational costs.
#' 
#' @note This function requires the \code{INLA} package, which is not a CRAN
#' package.  See \url{https://www.r-inla.org/download-install} for easy installation instructions.
#' @author David Bolin \email{davidbolin@@gmail.com} and Finn Lindgren \email{finn.lindgren@@gmail.com}
#' @references Bolin, D. and Lindgren, F. (2015) \emph{Excursion and contour uncertainty regions for latent Gaussian models}, JRSS-series B, vol 77, no 1, pp 85-106.
#' 
#' Bolin, D. and Lindgren, F. (2018), \emph{Calculating Probabilistic Excursion Sets and Related Quantities Using excursions}, Journal of Statistical Software, vol 86, no 1, pp 1-20.
#' @seealso \code{\link{excursions}}, \code{\link{excursions.mc}}
#'
#' @examples
#' ## In this example, we calculate the excursion function
#' ## for a partially observed AR process.
#' \dontrun{
#' if (require.nowarnings("INLA")) {
#' ## Sample the process:
#' rho = 0.9
#' tau = 15
#' tau.e = 1
#' n = 100
#' x = 1:n
#' mu = 10*((x<n/2)*(x-n/2) + (x>=n/2)*(n/2-x)+n/4)/n
#' Q = tau*sparseMatrix(i=c(1:n, 2:n), j=c(1:n, 1:(n-1)),
#'                      x=c(1,rep(1+rho^2, n-2),1, rep(-rho, n-1)),
#'                      dims=c(n, n), symmetric=TRUE)
#' X = mu+solve(chol(Q), rnorm(n))
#'
#' ## measure the sampled process at n.obs random locations
#' ## under Gaussian measurement noise.
#' n.obs = 50
#' obs.loc = sample(1:n,n.obs)
#' A = sparseMatrix(i=1:n.obs, j=obs.loc, x=rep(1, n.obs), dims=c(n.obs, n))
#' Y = as.vector(A %*% X + rnorm(n.obs)/sqrt(tau.e))
#'
#' ## Estimate the parameters using INLA
#' ef = list(c(list(ar=x),list(cov=mu)))
#' s.obs = inla.stack(data=list(y=Y), A=list(A), effects=ef, tag="obs")
#' s.pre = inla.stack(data=list(y=NA), A=list(1), effects=ef,tag="pred")
#' stack = inla.stack(s.obs,s.pre)
#' formula = y ~ -1 + cov + f(ar,model="ar1")
#' result = inla(formula=formula, family="normal", data = inla.stack.data(stack),
#'               control.predictor=list(A=inla.stack.A(stack),compute=TRUE),
#'               control.compute = list(config = TRUE,
#'                                      return.marginals.predictor = TRUE))
#'
#' ## calculate the level 0 positive excursion function
#' res.qc = excursions.inla(result, stack = stack, tag = 'pred', alpha=0.99, u=0,
#'                          method='QC', type='>', max.threads=2)
#' ## plot the excursion function and marginal probabilities
#' plot(res.qc$rho,type="l",
#'      main="marginal probabilities (black) and excursion function (red)")
#' lines(res.qc$F,col=2)
#' }}

excursions.inla <- function(result.inla,
                            stack,
                            name=NULL,
                            tag=NULL,
                            ind=NULL,
                            method,
                            alpha=1,
                            F.limit,
                            u,
                            u.link = FALSE,
                            type,
                            n.iter=10000,
                            verbose=0,
                            max.threads=0,
                            seed=NULL)
{
  if (!requireNamespace("INLA", quietly=TRUE))
    stop('This function requires the INLA package (see www.r-inla.org/download)')
  if(missing(result.inla))
	  stop('Must supply INLA result object')
  if(missing(method)){
	  cat('No method selected, using QC\n')
	  method = 'QC'
  }
  if(missing(u))
	  stop('Must specify level u')
  if(missing(type))
	  stop('Must specify type of excursion')

  if(result.inla$.args$control.compute$config==FALSE)
	  stop('INLA result must be calculated using control.compute$config=TRUE')

 if(missing(F.limit)) {
    F.limit = alpha
  } else {
    F.limit = max(alpha,F.limit)
  }

  n = length(result.inla$misc$configs$config[[1]]$mean)

  #Get indices for the component of interest in the configs
  ind.stack <- inla.output.indices(result.inla, name=name, stack=stack, tag=tag)
  n.out = length(ind.stack)
  #Index vector for the nodes in the component of interest
  ind.int <- seq_len(n.out)

  #ind is assumed to contain indices within the component of interest
  if(!missing(ind) && !is.null(ind)){
    ind <- private.as.vector(ind)
    ind.int <- ind.int[ind]
    ind.stack <- ind.stack[ind]
  }
  ind = ind.stack

  # If u.link is TRUE, the limit is given in linear scale
  # then transform to the scale of the linear predictor
  u.t = rho = rep(0,n)
  if(u.link == TRUE){
    links = result.inla$misc$linkfunctions$names[
                                            result.inla$misc$linkfunctions$link]
    u.tmp = sapply(ind, function(i) private.link.function(u,links[i]))
    u.t[ind] = u.tmp
  } else {
    u.t = u
  }

  #Calculate marginal probabilities
  #If stack and tag are provided, we are interested in the predictor
  #Otherwise, we are interested in some of the effects

  if(verbose)
	  cat('Calculating marginal probabilities\n')

  #Are we interested in a random effect?
  random.effect = FALSE
  if(!missing(name)&&!is.null(name)&&name!="APredictor"&&name!="Predictor"){
    random.effect = TRUE
    if(is.null(result.inla$marginals.random))
    stop('INLA result must be calculated using return.marginals.random=TRUE if excursion sets to be calculated for a random effect of the model')
  }

  if(!random.effect && is.null(result.inla$marginals.linear.predictor))
    stop('INLA result must be calculated using return.marginals.predictor=TRUE if excursion sets are to be calculated for the linear predictor.')

  if(random.effect) {
    rho.ind <- sapply(1:length(ind), function(j) inla.get.marginal(ind.int[j],
                           u=u,result = result.inla,
                           effect.name=name, u.link = u.link, type = type))
  } else {
    rho.ind <- sapply(1:length(ind), function(j) inla.get.marginal(ind[j],
                           u=u,result = result.inla, u.link = u.link, type = type))
  }
  rho[ind] = rho.ind

  n.theta = result.inla$misc$configs$nconfig
  for(i in 1:n.theta){
    config = private.get.config(result.inla,i)
    if(config$lp == 0)
      break
  }

  if(verbose)
    cat('Calculating excursion function using the ', method, ' method\n')

  if(method == 'EB' || method == 'QC' ) {
	  res = excursions(alpha=alpha, u=0, mu=config$mu-u.t, Q=config$Q,
	                   type=type, method=method, F.limit = F.limit,
	                   vars=config$vars, rho=rho,
		  	  		       ind=ind, n.iter=n.iter, max.threads=max.threads,seed=seed)
  	F = res$F[ind]
  } else if (method =='NI' || method == 'NIQC') {
    qc = 'QC'
  	if(method == 'NI')
	  	qc = 'EB'

  	res = lw = NULL

	  for(i in 1:n.theta){
		  if(verbose)
			  cat('Configuration ', i, ' of ', n.theta, '\n')

  		conf.i = private.get.config(result.inla,i)
	  	lw[i] = conf.i$lp
		  res[[i]] = excursions(alpha=alpha,u=0,mu=conf.i$mu-u.t,Q=conf.i$Q,
		                        type=type,method=qc,rho=rho,vars=conf.i$vars,
							              ind=ind,n.iter=n.iter,
							              F.limit = F.limit,
							              max.threads=max.threads,seed=seed)
  	}

  	w = exp(lw)/sum(exp(lw))

  	F = w[1]*res[[1]]$F[ind]
	  for(i in 2:n.theta){
 		  F = F + w[i]*res[[i]]$F[ind]
	  }
  } else if (method == 'iNIQC') {

  	pfam.i = rep(-0.1,n)
  	pfam.i[ind] = rho.ind
  	reo = sort(pfam.i,index.return=TRUE)$ix
  	pfam.i[!ind] = 0

	  res = lw =  NULL

	  for(i in 1:n.theta){
	  	if(verbose)
	  		cat('Configuration ', i, ' of ', n.theta, '\n')

  		conf.i = private.get.config(result.inla,i)
  		lw[i] = conf.i$lp
  		r.i <- INLA::inla(result.inla$.args$formula,
                        family = result.inla$.args$family,
                        data=result.inla$.args$data,
                        control.compute = list(config = TRUE,
                                               return.marginals.predictor = TRUE),
                        control.predictor=result.inla$.args$control.predictor,
                        control.mode = list(theta=
                         as.vector(result.inla$misc$configs$config[[i]]$theta),
                                  fixed=TRUE),
  		                  num.threads = "1:1")
  		# TODO: May refine the num.threads argument above to make it configurable, but
  		# since the inla() call is only constructing the model and optimising over the
  		# latent field once, for fixed hyperparameters, single threads is probably ok.

      if(random.effect) {
        p1.i <- sapply(1:length(ind), function(j) inla.get.marginal(ind.int[j],
                           u=u,result = r.i, effect.name=name,
                           u.link = u.link, type = type))

      } else {
        p1.i <- sapply(1:length(ind), function(j) inla.get.marginal(ind[j],
                           u=u,result = result.inla,
                           u.link = u.link, type = type))
      }
	  	pfam.i[ind] = p1.i

	  	res[[i]] = excursions(alpha=alpha,u=0,mu=conf.i$mu-u.t,
		  					            Q=conf.i$Q, type=type, method='QC',
	             						  rho=pfam.i,vars=conf.i$vars,
		  				          	  max.size=length(ind),reo=reo,
		  				          	  F.limit = F.limit,
			  				            n.iter=n.iter,max.threads=max.threads,seed=seed)
	  }

  	w = exp(lw)/sum(exp(lw))
  	F = w[1]*res[[1]]$F[ind]
  	for(i in 2:n.theta){
  		F = F + w[i]*res[[i]]$F[ind]
  	}
  } else {
  	stop('Method must be one of EB, QC, NI, NIQC, iNIQC')
  }

  if(type == "="){
    rho.ind = 1-pmax(rho.ind,1-rho.ind)
  }
  if(type == "!="){
    rho.ind = pmax(rho.ind,1-rho.ind)
  }

  F.out = mu.out = rho.out = vars.out = E.out = M.out = G.out = rep(NA,n.out)

  F.out[ind.int] = F
  vars.out[ind.int] = config$vars[ind]
  mu.out[ind.int] = config$mu[ind]
  rho.out[ind.int] = rho.ind


  G = rep(0,length(config$mu[ind]))
  if(type == "<") {
    G[config$mu[ind] > u.t[ind]] = 1
  } else {
    G[config$mu[ind] >= u.t[ind]] = 1
  }
  G.out[ind.int] = G

  E = rep(0,length(F))
  E[F>1-alpha] = 1
  E.out[ind.int] = E

  M = rep(-1,length(F))
  if (type=="<") {
    M[E==1] = 0
  } else if (type == ">") {
    M[E==1] = 1
  } else if (type == "!=" || type == "=") {
    M[E==1 & config$mu[ind]>u] = 1
    M[E==1 & config$mu[ind]<u] = 0
  }

  M.out[ind.int] = M

  output <- list(E = E.out,
                 F = F.out,
                 G = G.out,
                 M = M.out,
                 mean = mu.out,
                 vars = vars.out,
                 rho=rho.out,
                 meta=list(calculation="excursions",
                           type=type,
                           level=u,
                           level.link=u.link,
                           alpha=alpha,
                           F.limit = F.limit,
                           n.iter=n.iter,
                           method=method,
                           ind=ind.int,
                           call = match.call()))
  class(output) <- "excurobj"
  output
}
