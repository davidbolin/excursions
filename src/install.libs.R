mac = length(grep("mac", .Platform$pkgType)) > 0

if(mac){
	#change to libgfortran in Rhome:
	libs <- system(paste(c("otool -L ", "excursions"),collapse=""),intern=TRUE) 
	f.lib <- libs[grepl("libgfortran",libs)]
	f.lib2 = substr(f.lib,gregexpr("/",f.lib)[[1]][1],
					gregexpr(" ",f.lib)[[1]][1]-1)
	Rf.path <- paste(c(R.home(),"/lib/"),collapse="")
	Rlibs <- system(paste(c("ls ", Rf.path),collapse=""),intern=TRUE)
	Rf.lib <- Rlibs[grepl("libgfortran",Rlibs)]

	if(grepl("usr/local/lib/",f.lib) && grepl(Rf.lib,f.lib)){ 
		system(paste(c("install_name_tool -change ", "'",f.lib2, "' '",
			Rf.path, Rf.lib,"' excursions"),collapse=""))
		system(paste(c("install_name_tool -change ", "'",f.lib2, "' '",
			Rf.path, Rf.lib,"' gaussint"),collapse=""))
	}
	
	dest <- file.path(R_PACKAGE_DIR,"bin/mac")
	source.file1 <- "excursions"
	source.file2 <- "gaussint"	
	
} else if((.Platform$OS.type == "unix") && !mac){

	dest <- file.path(R_PACKAGE_DIR,"bin/linux")
	source.file1 <- "excursions"
	source.file2 <- "gaussint"
	
} else if(.Platform$OS.type == "windows") {
	#install to R_ARCH subdirectory:
	#if (nzchar(R_ARCH)){
	#	libarch <- paste('bin/windows', R_ARCH, sep='') 
	#} else {
	#	libarch <- "bin/windows"
	#}
	if (R.Version()$arch == "i386"){
		libarch <- "bin/windows/i386" 
	} else {
		libarch <- "bin/windows/x64"
	}
	dest <- file.path(R_PACKAGE_DIR,libarch)

	dest <- file.path(R_PACKAGE_DIR,libarch)
	source.file1 <- "excursions.exe"
	source.file2 <- "gaussint.exe"
	
}else {
	stop("OS not supported")
}
dir.create(dest,recursive=TRUE,showWarnings=FALSE)

file.copy(source.file1,dest,overwrite=TRUE) 
file.copy(source.file2,dest,overwrite=TRUE) 





