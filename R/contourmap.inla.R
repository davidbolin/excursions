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


contourmap.inla <- function(result.inla, stack,name=NULL,tag=NULL, method="EB",
                            n.levels,ind,levels, type='standard', measure,
                            contour.map.function = FALSE, use.marginals=TRUE,
                            alpha=1,n.iter=10000, verbose=FALSE,max.threads=0)
{
  if (!require("INLA"))
    stop('This function requires the INLA package (see www.r-inla.org/download)')
  if(missing(result.inla))
    stop('Must supply INLA result object')

  if(method != "EB")
    stop("only EB method implemented so far")

  if(verbose) cat("calculating contour map using the EB method")

  ind.stack <- inla.output.indices(result.inla, name=name, stack=stack, tag=tag)
  ind.int <- seq_len(length(ind.stack))
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
  return(contourmap(mu=config$mu,Q = config$Q, n.levels = n.levels,
                    ind = ind,levels = levels, type=type, measure = measure,
                    contour.map.function = contour.map.function,
                    use.marginals=use.marginals, alpha=alpha, n.iter=n.iter,
                    verbose=verbose,max.threads=max.threads))
}