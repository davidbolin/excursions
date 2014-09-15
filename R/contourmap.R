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


contourmap <- function(mu,Q,vars,n.levels,ind,levels,
								       type='standard', measure,
								       contour.map.function = FALSE,
								       use.marginals=TRUE,
								       alpha=1,n.iter=10000,
								       verbose=FALSE,max.threads=0)
{
	if(type == 'standard')
	{
	  if(verbose) cat('Creating contour map\n')
		lp <- excursions.levelplot(mu=mu,n.levels = n.levels,ind = ind,
		                           levels = levels,equal.area=FALSE)
	}
	else if(type == 'equalarea')
	{
		if(verbose) cat('Creating equal area contour map\n')
		lp <- excursions.levelplot(mu = mu,n.levels = n.levels,ind = ind,
		                           levels = levels,equal.area=TRUE)
	}
	else if(type == 'P0-optimal' || type == 'P1-optimal' || type == 'P2-optimal')
  {
		if(!missing(levels)){
			warning('Not using supplied levels for optimal contour map\n')
			if(!missing(n.levels)){
				if(n.levels != length(levels)){
					warning('n.levels != length(levels), using n.levels\n')
				}
			} else {
				n.levels = length(levels)
			}
		}
		if(missing(vars) && missing(Q)){
			stop('Variances must be supplied when creating optimal contour map')
    } else if(missing(vars)) {
      L = chol(Q)
      vars = excursions.variances(L)
    }
    if(use.marginals == TRUE){
			if(missing(Q))
			  stop('The precision matrix must be supplied unless marginals are used')
		}

		if(type == 'P0-optimal'){
			if(verbose) cat('Creating P0-optimal contour map\n')
			opt.measure = 0
		}else if(type == 'P1-optimal'){
			if(verbose) cat('Creating P1-optimal contour map\n')
			opt.measure = 1
		} else if (type == 'P2-optimal'){
			if(verbose) cat('Creating P2-optimal contour map\n')
			opt.measure = 2
		} else {
			stop('only P0, P1, and P2-measures supported')
		}

		lp <- excursions.opt.levelplot(mu = mu,vars = vars,Q = Q,
		                               n.levels = n.levels, measure = opt.measure,
		                               use.marginals = use.marginals,ind = ind)
	} else {
		stop('type not supported (use standard, equalarea, or P0/P1/P2-optimal)')
	}

	F.calculated = FALSE
	if(!missing(measure) && !is.null(measure)){
		if(missing(Q))
			stop('precision matrix must be supplied if measure should be calculated')

		for( i in 1:length(measure)){
      if(measure[i]==1) {
        if(verbose) cat('Calculating P1-measure\n')
        lp$P1 <- Pmeasure(lp=lp,mu=mu,Q=Q,ind=ind,type=measure[i])
      } else if(measure[i] == 2) {
        if(verbose) cat('Calculating P2-measure\n')
        lp$P2 <- Pmeasure(lp=lp,mu=mu,Q=Q,ind=ind,type=measure[i])
	 	  } else if (measure[i] == 0) {
	 	    if(verbose) cat('Calculating P0-measure and contour map function\n')

	 	    p <- contourfunction(lp=lp, mu=mu,Q=Q ,vars=vars, ind = ind,
	 	                         alpha=alpha, n.iter=n.iter,max.threads=max.threads)

	 	    lp$P0 = mean(p$F)
	 	    lp$F = p$F
	 	    F.calculated = TRUE
	 	  } else {
	 	    cat('Measure not supported (Only 0,1, and 2 supported)\n')
	 	  }
	 	}
	}
	if(contour.map.function == TRUE && F.calculated == FALSE){
    if(verbose) cat('Calculating contour map function\n')
    p <- contourfunction(lp=lp, mu=mu,Q=Q ,vars=vars, ind = ind,
	 	                     alpha=alpha, n.iter=n.iter, max.threads=max.threads)
	 	lp$P0 = mean(p$F)
	 	lp$F = p$F
	}
	lp$meta <- list(calculation="contourmap",
                        alpha=alpha,
                        n.iter=n.iter,
                        levels=lp$u,
                        alpha=alpha,
                        type="!=")
        class(lp) <- "excurobj"
	return(lp)
}
