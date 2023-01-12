## Find the indices into inla output config structures corresponding
## to a specific predictor, effect name, or inla.stack tag.
##
## result : an inla object
inla.output.indices = function(result, name=NULL, stack=NULL, tag=NULL, ...)
{
  
    if (!is.null(result$misc$configs$.preopt) && result$misc$configs$.preopt) {
      inla.experimental <- TRUE
    } else {
      inla.experimental <- FALSE  
    }
    
    if (!is.null(name) && !is.null(tag)) {
        stop("At most one of 'name' and 'tag' may be non-null.")
    }
    if (!is.null(tag)) {
        if (is.null(stack)) {
            stop("'tag' specified but 'stack' is NULL.")
        }
        tags <- names(stack$data$index)
        if (!(tag %in% tags)) {
            stop("'tag' not found in 'stack'.")
        }
    } else if (is.null(name) || (!is.null(name) && (name == ""))) {
        if ("APredictor" %in% result$misc$configs$contents$tag) {
            name <- "APredictor"
        } else {
            name <- "Predictor"
        }
    }
    
    if (inla.experimental &&
        !is.null(name) &&
        (name == "APredictor" || name == "Predictor") ){
        stop("INLA was run in experimental mode, so you can only 
           compute excursion sets for model components. This may be improved in a future version.")
    }
    
    ## Find variables
    if (!is.null(name)) {
        if (!(name %in% result$misc$configs$contents$tag)) {
            stop("'name' not found in result.")
        }
        ct <- result$misc$configs$contents
        
        #Shift indices for experimental mode
        if (inla.experimental){ 
            for(nm in c("APredictor", "Predictor")) {
                if (ct$tag[1] == nm) {
                    ct$tag <- ct$tag[-1]
                    ct$start <- ct$start[-1] - ct$start[2] + 1
                    ct$length <- ct$length[-1]
                }
            }
        }
        
        nameindex <- which(ct$tag == name)
        index <- (ct$start[nameindex] - 1L + seq_len(ct$length[nameindex]))  
        
    } else { ## Have tag
        if (inla.experimental) {
            stop("INLA was run in experimental mode, so you can only 
           compute excursion sets for model components. This may be improved in a future version.")
        }
        # FL 2023-01-12: the internal latent mean is in improved.mean for each
        # configuration, and the A-matrices are available:
        # result$misc$configs$pA %*% result$misc$configs$A
        # Together with the Qinv info in each configuration, this should cover the needed
        # information.
        index <- stack$data$index[[tag]]
    }

    index
}

private.simconf.link <- function(res,links,trans=TRUE)
{
  if(trans){
    n = length(res$a)
    res$a.marginal = sapply(1:n, function(i) private.link.function(
                                      res$a.marginal[i],links[i],inv=TRUE))
    res$b.marginal = sapply(1:n, function(i) private.link.function(
                                      res$b.marginal[i],links[i],inv=TRUE))
    res$a = sapply(1:n, function(i) private.link.function(
                                      res$a[i],links[i],inv=TRUE))
    res$b = sapply(1:n, function(i) private.link.function(
                                      res$b[i],links[i],inv=TRUE))
  }
  return(res)
}

private.link.function <- function(x, link, inv=FALSE)
{
  if (is.na(link)) {
    link = "identity"
  }
  return(do.call(paste("inla.link.", link, sep=""),list(x=x, inv=inv)))
}

private.get.config <- function(result,i)
{
  mu=result$misc$configs$config[[i]]$mean
  Q=forceSymmetric(result$misc$configs$config[[i]]$Q)
  vars = diag(result$misc$configs$config[[i]]$Qinv)
  m <- max(unlist(lapply(result$misc$configs$config,
                         function(x) x$log.posterior)))
  lp = result$misc$configs$config[[i]]$log.posterior -m
  return(list(mu=mu,Q=Q,vars=vars,lp=lp))
}

## Calculate the marginal probability for X_i>u or X_i<u.
## Note that the index 'i' refers to a location in the linear
## predictor if predictor==TRUE, whereas it refers to a location
## in the random effect vector otherwise.
inla.get.marginal <- function(i, u,result,effect.name=NULL, u.link, type)
{
  if(is.null(effect.name) && u.link == TRUE){
      marg.p = result$marginals.fitted.values[[i]]
  } else if (is.null(effect.name)){
      # Calculate marginals using linear predictor
      marg.p = result$marginals.linear.predictor[[i]]
  } else {
      # Calculate marginals using a random effect
      marg.p = result$marginals.random[[effect.name]][[i]]
  }

  if(type=='<'){
	  return(INLA::inla.pmarginal(u,marg.p))
  } else {
	  return(1-INLA::inla.pmarginal(u,marg.p))
  }
}

## Calculate the marginal probability for a<X_i<b.
## The function returns c(P(X<a),P(X<b))
inla.get.marginal.int <- function(i, a,b,result,effect.name=NULL)
{
  if (is.null(effect.name)){
    # Calculate marginals using linear predictor
    marg.p = result$marginals.linear.predictor[[i]]
  } else {
    # Calculate marginals using a random effect
    marg.p = result$marginals.random[[effect.name]][[i]]
  }
  return(c(INLA::inla.pmarginal(a,marg.p),INLA::inla.pmarginal(b,marg.p)))

}