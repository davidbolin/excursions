# Description #
The package contains functions that compute probabilistic excursion sets, contour credibility regions, and simultaneous confidence bands for latent Gaussian random processes and fields. 

The theory of the methods used in the package are described in [Bolin and Lindgren (2014)](http://onlinelibrary.wiley.com/doi/10.1111/rssb.12055/abstract). A manual for the package can be found [here](http://www.math.chalmers.se/~bodavid/software/excursions/excursions_manual.pdf) and the code used in the manual can be downloaded [here](http://www.math.chalmers.se/~bodavid/software/excursions/code.zip).

# Version #
The current stable version of the package is 1.1, which is also on CRAN. 

The development version, which is called excursionsdevel, can be installed directly from R using the command 
```
#!r

devtools::install_bitbucket("excursions","davidbolin",ref="devel")
```
