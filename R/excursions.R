## excursions.R
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

excursions <- function(alpha, u, mu, Q, type, n.iter=10000, Q.chol,
                       vars, rho, reo, method='EB', ind, max.size,
                       verbose=0, keep=FALSE, max.threads=0) {

# Paths to temp directoy
tmppath = private.make.tempdir()

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Extract input data and setup parameters %
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(method=='QC'){
	qc = TRUE
} else if(method == 'EB'){
	qc = FALSE
} else {
	stop('only EB and QC methods are supported.')
}
if(missing(alpha))
	stop('Must specify error probability')

if(missing(u))
	stop('Must specify level')

if(missing(mu))
	stop('Must specify mean value')

if(missing(Q) && missing(Q.chol))
	stop('Must specify a precision matrix or its Cholesky factor')

if(missing(type))
	stop('Must specify type of excursion set')

if(qc && missing(rho))
	stop('rho must be provided if QC is used.')

if(!missing(ind) && !missing(reo))
	stop('Either provide a reordering using the reo argument or provied a set of nodes using the ind argument, both cannot be provided')


if (missing(rho)) {
	rho_p = 0
} else {
	if(qc) {
		rho_p = 2
	} else {
		rho_p = 1
	}
}

if (missing(reo)){
	reo_p = 0
} else {
	reo_p = 1
}

if (missing(vars)) {
	var_p = 0
} else {
	var_p = 1
}

if (missing(max.size)){
	m.size = length(mu)
} else {
	m.size = max.size
}
if (missing(ind)) {
	ind_p = 0
} else {
	ind_p = 1
	indices = rep(0,length(mu))
	indices[ind] = 1
	if(missing(max.size)){
		m.size = length(ind)
	} else {
		m.size = min(length(ind),m.size)
	}
}



if (!missing(Q.chol) && !is.null(Q.chol)) {
    ## make the representation unique (i,j,v)
    Q = private.as.dgTMatrix(Q.chol)
    is.chol = TRUE
} else {
    ## make the representation unique (i,j,v)
    Q = private.as.dgTMatrix(Q)
    is.chol = FALSE
}



if(type == '>'){
	tp = 0
} else if(type == '<') {
	tp = 1
} else if(type == '=' || type == '!='){
	tp = 2
} else {
	stop('Type must be one of <, >, =, or !=.')
}

#%%%%%%%%%%%%%%%%%%%%%
# Write data to disk %
#%%%%%%%%%%%%%%%%%%%%%

fid <- file(file.path(tmppath,"initdata.bin"), "wb")
writeBin(as.double(alpha),fid)
writeBin(as.double(u), fid)
writeBin(as.integer(n.iter), fid)
writeBin(as.integer(tp), fid)
writeBin(as.integer(is.chol),fid)
writeBin(as.integer(rho_p), fid)
writeBin(as.integer(reo_p), fid)
writeBin(as.integer(var_p), fid)
writeBin(as.integer(ind_p), fid)
writeBin(as.integer(m.size), fid)
writeBin(as.integer(verbose), fid)
writeBin(as.integer(max.threads), fid)
close(fid);

fid <- file(file.path(tmppath,"mu.bin"), "wb")
writeBin(mu,fid);
close(fid);

## write upper triangular and and diagonal part of Q to file:
uppertri = (Q@i <= Q@j)
i = Q@i[uppertri]
j = Q@j[uppertri]
v = Q@x[uppertri]

## dgTMatrix internals are sorted by j; sort by j instead:
is = sort(i,index.return=TRUE)
i = is$x
j = j[is$ix]
v = v[is$ix]
fid <- file(file.path(tmppath,"precI.bin"), "wb")
writeBin(as.integer(c(matrix(c(i,j),nrow = 2,byrow=T))), fid)
close(fid)
fid <- file(file.path(tmppath,"precV.bin"), "wb")
writeBin(as.double(v), fid)
close(fid)

if (rho_p >0){
	fid <- file(file.path(tmppath,"rho.bin"), "wb")
	writeBin(as.double(rho), fid)
	close(fid)
}
if (var_p == 1){
	fid <- file(file.path(tmppath,"vars.bin"), "wb")
	writeBin(as.double(vars), fid)
	close(fid)
}
if (reo_p == 1){
	fid <- file(file.path(tmppath,"reo.bin"), "wb")
	writeBin(as.integer(reo-1), fid)
	close(fid)
}

if (ind_p == 1){
	fid <- file(file.path(tmppath,"ind.bin"), "wb")
	writeBin(as.integer(indices), fid)
	close(fid)
}

#%%%%%%%%%%%%%%%%%%%%%%
# Do the calculations %
#%%%%%%%%%%%%%%%%%%%%%%
private.excursions.call(gsub("/*$", "/", tmppath),"excursions")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Read results and save data %
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = length(mu)

fid <- file(file.path(tmppath,"results_p.bin"), "rb")
Pv = readBin(fid,double(),n)
close(fid)

fid <- file(file.path(tmppath,"results_e.bin"), "rb")
Ev = readBin(fid,double(),n)
close(fid)

fid <- file(file.path(tmppath,"results_rho.bin"), "rb")
rho = readBin(fid,double(),n)
close(fid)

fid <- file(file.path(tmppath,"results_reo.bin"), "rb")
ind = readBin(fid,integer(),n)
close(fid)

fid <- file(file.path(tmppath,"vars.bin"), "rb")
vars = readBin(fid,double(),n)
close(fid)

if (!keep) {
    ## Delete temporary files
    unlink(tmppath, recursive=TRUE)
}

ind = ind + 1

ii = which(Pv[1:n] > 0)
if (length(ii) == 0) i=n+1 else i=min(ii)

F = Fe  = rep(0,n)
F[ind] = Pv
Fe[ind] = Ev

ireo = NULL
ireo[ind] = 1:n

if(type == '='){
	F=1-F
	D = rep(1,n)
	if(i<n+1) D[ind[i:n]] = 0
} else {
	D = rep(0,n)
	if(i<n+1) D[ind[i:n]] = 1
}

if (keep) {
    return(list(F=F, Fe=Fe, D=D, rho=rho, reo=ind, ireo=ireo, vars=vars,
                tmppath=file.path(tmppath, "")))
} else {
    return(list(F=F, Fe=Fe, D=D, rho=rho, reo=ind, ireo=ireo, vars=vars))
}
}
