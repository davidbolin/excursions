## utils.R 
##
##   Copyright (C) 2013 Håvard Rue, David Bolin, Finn Lindgren
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

.onLoad <- function(lib,pkg){
	#Set environment variable to find BLAS, LAPACK and MATRIX on windows
	if(.Platform$OS.type == "windows"){
		rhome = Sys.getenv("R_HOME")
		rbin = paste(c(rhome,"/bin/",.Platform$r_arch),collapse="")
		path = Sys.getenv("PATH")
		newpath = paste(c(path,rbin),collapse=";")
		Sys.setenv(PATH=newpath)
	}
}


private.link.function <-
    function(x, link, inv=FALSE)
{
    if (is.na(link)) {
        link = "identity"
    }
    return(do.call(paste("inla.link.", link, sep=""),list(x=x, inv=inv)))
}

private.excursions.call = function (path,func)
{
	if(private.os.type() == "mac"){
		binpath = file.path(system.file("bin",package="excursions"),"mac")
	} else if(private.os.type() == "linux"){
		binpath = file.path(system.file("bin",package="excursions"),"linux")
	} else if(private.os.type() == "windows"){
		if(R.Version()$arch == "i386"){
			binpath = file.path(system.file("bin",package="excursions"),
												"windows/i386")
		} else {
			binpath = file.path(system.file("bin",package="excursions"),
												"windows/x64")
		}
	} else {
		stop("OS type not supported")
	}
	system(paste(file.path(binpath,func), file.path(path,"")))
}

#from INLA
private.as.dgTMatrix = function (A, unique = TRUE)
{
    if (unique) {
        return(as(as(as(A, "CsparseMatrix"), "dgCMatrix"), "dgTMatrix"))
    }
    else {
        if (is(A, "dgTMatrix")) {
            return(A)
        }
        else {
            return(as(as(A, "TsparseMatrix"), "dgTMatrix"))
        }
    }
}

private.tempfile = function (pattern = "file", tmpdir = tempdir())
{
    return(gsub("\\\\", "/", tempfile(pattern, tmpdir)))
}
private.tempdir = function ()
{
    return(gsub("\\\\", "/", tempdir()))
}

private.make.tempdir = function (dir = NULL)
{
    if (is.null(dir)) {
        dir = private.tempfile("exctmp")
    }
    dir.create(dir)
    return(dir)
}

private.binpath = function ()
{
    return("../excursions_v2/")
}

private.os = function(type = c("linux", "mac", "windows", "else"))
{
    if (missing(type))
        stop("Type of OS is required.")
    type = match.arg(type)

    if (type == "windows")
        return (.Platform$OS.type == "windows")
    else if (type == "mac")
        return (length(grep("mac", .Platform$pkgType)) > 0)
    else if (type == "linux")
        return ((.Platform$OS.type == "unix") && !private.os("mac"))
    else if (type == "else")
        return (TRUE)
    else
        stop("This shouldn't happen.")
}
private.os.type = function()
{
    for (os in c("windows", "mac", "linux", "else"))
        if (private.os(os))
            return (os)
    stop("This shouldn't happen.")
}






