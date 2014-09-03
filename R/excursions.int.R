## excursions.int.R
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


excursions.int <- function(mu,Q,a,b,n.iter=10000,Q.chol,ind,verbose=0,alpha=1,max.threads=0){

# Paths to temp directoy
tmppath = private.make.tempdir()

# Extract input data and setup parameters
if(missing(mu))
	stop('Must specify mean value')

if(missing(Q) && missing(Q.chol))
	stop('Must specify a precision matrix or its Cholesky factor')

if(missing(a))
	stop('Must lower limit')

if(missing(b))
	stop('Must upper limit')

if(!missing(ind) && !missing(Q.chol))
	stop('Cannot provide both cholesky factor and indices.')

m.size = length(mu)

if (missing(ind) || is.null(ind)) {
	ind_p = 0
} else {
	ind_p = 1
	indices = rep(0,length(mu))
	indices[ind] = 1
	m.size = length(ind)
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

# Write data to disk
fid <- file(file.path(tmppath,"initdata.bin"), "wb")
writeBin(as.double(alpha),fid)
writeBin(as.integer(n.iter), fid)
writeBin(as.integer(is.chol),fid)
writeBin(as.integer(ind_p), fid)
writeBin(as.integer(m.size), fid)
writeBin(as.integer(verbose), fid)
writeBin(as.integer(max.threads), fid)
close(fid);

fid <- file(file.path(tmppath,"mu.bin"), "wb")
writeBin(mu,fid);
close(fid);

fid <- file(file.path(tmppath,"a.bin"), "wb")
writeBin(a,fid);
close(fid);

fid <- file(file.path(tmppath,"b.bin"), "wb")
writeBin(b,fid);
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

if (ind_p == 1){
	fid <- file(file.path(tmppath,"ind.bin"), "wb")
	writeBin(as.integer(indices), fid)
	close(fid)
}

# Do the calculations
private.excursions.call(gsub("/*$", "/", tmppath),"gaussint")

# Read results and save data
fid <- file(file.path(tmppath,"results.bin"), "rb")
res = readBin(fid,double(),2)
close(fid)

## Delete temporary files
unlink(tmppath, recursive=TRUE)

return(list(P=res[1],Pe=res[2]))
}
