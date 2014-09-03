## excursions.inla.R
##
##   Copyright (C) 2012, 2013 David Bolin, Finn Lindgren
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

excursions.inla <- function(result.inla, ind=NULL, method, alpha=1, u,
                            type, n.iter=10000, verbose=0, max.threads=0){
#Wrapper function that takes INLA output and calculates excursion function
if (!require("INLA")) {
    stop('This function requires the INLA package (see www.r-inla.org/download)')
}

#ind <- inla.stack.index(stack, 'pred')$data
if(missing(result.inla))
	stop('Must supply INLA result object')
if(missing(method)){
	cat('No method selected, using QC')
	method = 'QC'
}
if(missing(u))
	stop('Must specify level u')
if(missing(type))
	stop('Must specify type of excursion')

if(result.inla$.args$control.compute$config==FALSE)
	stop('INLA result must be calculated using control.compute$config=TRUE')

if(missing(ind)||is.null(ind))
	stop('You must provide indices of nodes in the joint distribution for which the excursion function should be calculated using the ind-argument')

#extract the link function
#tmp = INLA:::inla.models.section.likelihood()
#fam = result.inla$.args$family
#link = result.inla$.args$link
#if(is.null(link))
#	link = tmp$likelihood[[fam]]$link[2]

links = result.inla$misc$linkfunctions$names[
				result.inla$misc$linkfunctions$link]



#transform limit
#u.t = private.link.function(u,link)
u.t = rep(0,length(result.inla$misc$configs$config[[1]]$mean))
u.tmp = sapply(ind, function(i) private.link.function(u,links[i]))
u.t[ind] = u.tmp

#calculate parameteric family
pfam = rep(0,length(result.inla$misc$configs$config[[1]]$mean))
if(verbose){
	cat('Calculating parametric family based on marginals\n')
}

for (i in ind){
	if(links[i]=="identity" || is.na(links[i])){
	marg.p = result.inla$marginals.fitted.values[[i]]
	} else {
	marg.p =
            INLA::inla.tmarginal(function(x)
                                 private.link.function(x,links[i],inv=TRUE),
                                 result.inla$marginals.fitted.values[[i]])
	}
	if(type=='<'){
		pfam[i] = INLA::inla.pmarginal(u,marg.p)
	}else {
		pfam[i] = 1-INLA::inla.pmarginal(u,marg.p)
	}
}

p1 = pfam[ind]
#p1 <- sapply(ind, function(i) 1- INLA::inla.pmarginal(u,INLA::inla.tmarginal(function(x)
#			private.link.function(x,links[i],inv=TRUE),
#			result.inla$marginals.fitted.values[[i]])))
#pfam[ind] = p1

if(verbose){
	cat('Find the mode and calculate variances\n')
}
n.theta = result.inla$misc$configs$nconfig
for(i in 1:n.theta){
	if (result.inla$misc$configs$config[[i]]$log.posterior == 0){
		mu.p=result.inla$misc$configs$config[[i]]$mean
		Q.p=result.inla$misc$configs$config[[i]]$Q
		#vars = diag(INLA::inla.qinv(Q.p))
		vars = diag(result.inla$misc$configs$config[[i]]$Qinv)
	}
}

marginal.prob = p1
mean.field = mu.p[ind]

if(method == 'EB' || method == 'QC' ){
	if(method == 'EB'){
		if(verbose){
			cat('Calculating excursion function using the EB method\n')
		}
		qc = 'EB'
	} else {
		if(verbose){
			cat('Calculating excursion function using the QC method\n')
		}
		qc = 'QC'
	}
	res = excursions(alpha=alpha, u=0, mu=mu.p-u.t, Q=Q.p, type=type,
					method=qc, vars=vars, rho=pfam, ind=ind, n.iter=n.iter,
					max.threads=max.threads)
	F = res$F[ind]

} else if (method =='NI' || method == 'NIQC') {
	if(method == 'NI'){
		if(verbose){
			cat('Calculating excursion function using the NI method\n')
		}
		qc = 'EB'
	} else {
		if(verbose){
			cat('Calculating excursion function using the NIQC method\n')
		}
		qc = 'QC'
	}
	res = NULL
	lw = NULL
	for(i in 1:n.theta){
		if(verbose){
			cat('Configuration ')
			cat(i)
			cat(' of ')
			cat(n.theta)
			cat('\n')
		}
		mu.p=result.inla$misc$configs$config[[i]]$mean
		Q.p=result.inla$misc$configs$config[[i]]$Q
		#vars = diag(INLA::inla.qinv(Q.p))
		vars = diag(result.inla$misc$configs$config[[i]]$Qinv)
		lw[i] = result.inla$misc$configs$config[[i]]$log.posterior
		res[[i]] = excursions(alpha=alpha,u=0,mu=mu.p-u.t,Q=Q.p,type=type,
							  method=qc,rho=pfam,vars=vars,
							  ind=ind,n.iter=n.iter,max.threads=max.threads)
	}

	w = exp(lw)/sum(exp(lw))

	#calculate the excursion function
	F = w[1]*res[[1]]$F[ind]
	for(i in 2:n.theta){
 		F = F + w[i]*res[[i]]$F[ind]
	}

} else if (method == 'iNIQC') {
	if(verbose){
		cat('Calculating excursion function using the improved NIQC method\n')
	}
	pfam.i = rep(-0.1,length(result.inla$misc$configs$config[[1]]$mean))
	pfam.i[ind] = p1
	pfam.s = sort(pfam.i,index.return=TRUE)
	reo = pfam.s$ix
	ireo = NULL
	ireo[reo] = 1:length(reo)
	pfam.i = rep(0,length(result.inla$misc$configs$config[[1]]$mean))
	pfam.i[ind] = p1

	res = NULL
	lw = NULL
	for(i in 1:n.theta){
		if(verbose){
			cat('Configuration ')
			cat(i)
			cat(' of ')
			cat(n.theta)
			cat('\n')
		}
		mu.p=result.inla$misc$configs$config[[i]]$mean
		Q.p=forceSymmetric(result.inla$misc$configs$config[[i]]$Q)
		#vars = diag(INLA::inla.qinv(Q.p))
		vars = diag(result.inla$misc$configs$config[[i]]$Qinv)
		lw[i] = result.inla$misc$configs$config[[i]]$log.posterior
		r.i <- INLA::inla(result.inla$.args$formula,
                                  family = result.inla$.args$family,
                                  data=result.inla$.args$data,
                                  control.compute = list(config = TRUE),
                                  control.predictor=result.inla$.args$control.predictor,
                                  control.mode = list(theta=
                                  as.vector(result.inla$misc$configs$config[[i]]$theta),
                                  fixed=TRUE))

		p1.i <- sapply(ind, function(j) 1-INLA::inla.pmarginal(u,
					INLA::inla.tmarginal(function(x)
									private.link.function(x,links[i],inv=TRUE),
									r.i$marginals.fitted.values[[j]])))
		pfam.i[ind] = p1.i
		res[[i]] = excursions(alpha=alpha,u=0,mu=mu.p-u.t,
							Q=Q.p, type=type, method='QC',
							  rho=pfam.i, vars=vars,
							  max.size=length(ind),reo=reo,
							  n.iter=n.iter,max.threads=max.threads)
	}

	w = exp(lw)/sum(exp(lw))
	F = w[1]*res[[1]]$F[ind]
	for(i in 2:n.theta){
 		F = F + w[i]*res[[i]]$F[ind]
	}
} else {
	stop('Method must be one of EB, QC, NI, NIQC, iNIQC')
}

return(list(F=F,mean=mean.field,rho=marginal.prob))
}

