# Description #
The package contains functions that compute probabilistic excursion sets, contour credibility regions, and simultaneous confidence bands for latent Gaussian random processes and fields. 

The theory of the methods used in the package are described in [Bolin and Lindgren (2014)](http://onlinelibrary.wiley.com/doi/10.1111/rssb.12055/abstract). A manual for the package can be found [here](http://www.math.chalmers.se/~bodavid/software/excursions/excursions_manual.pdf) and the code used in the manual can be downloaded [here](http://www.math.chalmers.se/~bodavid/software/excursions/code.zip).

# Version #
The current stable version of the package is 1.1, which is also on CRAN. 

The development version can be installed directly from R using the command 
```
#!r

devtools::install_bitbucket("excursions","davidbolin",ref="devel")
```
The development version of the package is called excursionsdevel, so it is loaded by running

```
#!r
library(excursionsdevel)
```

If you want to install the development version on Windows, you first need to install Rtools and update the Windows PATH environment variable. This can be done for the current R session only 
```
#!r
rtools = "C:\\Rtools\\bin"
gcc = "C:\\Rtools\\gcc-4.6.3\\bin"
Sys.setenv(PATH=paste(c(gcc,rtools, Sys.getenv("PATH")),collapse=";"))
```
where the variables rtools and gcc need to be changed if Rtools is not installed directly on C:. 