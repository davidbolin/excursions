# Development version

* Fix `excursions.inla()`, `contourmap.inla()` and `simconf.inla()` for
  linear predictors in compact INLA mode: offsets are now included in the
  predictor mean, and the linear predictor is linked to the latent field
  with the INLA default precision `exp(15)` instead of `1e9`, which made
  the joint precision numerically singular (`D[i,i] is negative`) for
  larger models.
* Fix the marginal variances of linear predictors in compact INLA mode in
  `excursions.inla()`, `contourmap.inla()` and `simconf.inla()`. They were
  computed from the sparse partial inverse stored by INLA, which misses
  covariances between latent components, and are now computed from the joint
  precision. The wrong variances made `method = "QC"` inconsistent (excursion
  function larger than the marginal probabilities) and gave suboptimal
  orderings for the other methods.
* Fix `continuous()` for excursion and contour map objects computed with an
  `ind` argument that selects a strict subset of the geometry nodes, which
  failed with a `non-conformable arguments` error.
* Faster marginal variances in `excursions.variances()`, and therefore in
  most other functions, using a scatter-based Takahashi recursion on the
  pattern of the Cholesky factor. The argument `max.threads` of 
  `excursions.variances()` is no longer used.
* Faster sequential importance sampling in `gaussint()`, `excursions()`,
  `contourmap()` and `simconf()`. With several threads, all of the work is
  now done in parallel, and the results for a given seed only depend on the
  number of threads. Results for a single thread are unchanged.
* Faster `tricontour()`, `tricontourmap()` and `continuous()`, which were
  quadratic in the number of contour segments, and faster `simconf.mc()`,
  `excursions.mc()`, `contourmap.mc()`, `simconf.mixture()` and
  `simconf.inla()` with `method = "NI"`.
* Much faster reordering of the nodes in `excursions()`, `contourmap()` and
  the INLA interfaces when many nodes have marginal probabilities above the
  threshold, by merging single node constraint sets before calling CAMD.
  With `excursions.inla()`, the joint variances of the linear predictor are
  now only computed for the configurations that are used. 
* Fix the reordering of the nodes above the threshold, where the two nodes
  with the largest marginal probabilities were in the same constraint set
  and could be integrated in the wrong order. They are now integrated in the
  order of their probabilities like the other nodes.
* The quantiles of Gaussian mixtures in `simconf.mixture()` and
  `simconf.inla()` are now computed to full precision.
* Fix `simconf()`, which returned `a.marginal` and `b.marginal` swapped, and
  ignored `n.iter`.
* Fix `gaussint()` and `simconf()`, which ignored `ind` when it was given
  as integer indices. This also affects `simconf.inla()` with `ind`.
* Fix `excursions()` and `contourmap()` with `Q.chol`, where the Cholesky
  factor was used as a precision matrix when the nodes were reordered.
* Fix `excursions.inla()` with `method = "iNIQC"`, which failed for models
  with fixed hyperparameters since the refits did not use all the arguments
  of the original fit, and which used the marginals of the original fit
  instead of the refits for the linear predictor.
* Fix `simconf.mixture()`: `mix.samp = FALSE` failed without `ind`, and
  permuted computed variances twice when `ind` was given; `seed` is now used
  for `mix.samp = TRUE`; the returned `mean` and `vars` are now the mean and
  variances of the mixture; and a single mixture component is supported.
* Fix a wrong quantile in `simconf.inla()` with `method = "NI"` when the
  limits of the mixture quantiles had to be extended.
* `simconf.mc()` no longer prints the estimated level.
* `excursions.variances()` accepts supernodal Cholesky factors.

# excursions 2.5.11
 
* Update link to external inla vignette in documentation

# excursions 2.5.10
 
* Corrections to various links in documentation 

# excursions 2.5.9

* Add -DNPRINT to CAMD build, and remove fflush(stdout) use, to avoid fprint
  and stdout usage
* Fix issue with identity matrices, ensuring they are not treated as all-zero
* Remove obsolete and unneded includes of Rdefines.h and R_ext/PrtUtil.h

# excursions 2.5.8

* Minor updates to C code to avoid warning with gcc13 and clang17 compilers

# excursions 2.5.7

* Update to make sure that max.threads properly limits the number of threads
* Minor update to the documentation

# excursions 2.5.6

* Remove rgdal suggest
* Update some examples to use fmesher instead of INLA

# excursions 2.5.5

* Add R compiler configuration extraction to configure.ac
* Add support for results computed with the new default `compact` mode in INLA

# excursions 2.5.4

* Remove local gsl code, and rely on SystemRequirements instead

# excursions 2.5.3

* Add support for INLA results computed with the `experimental` option.
* Avoid deprecated Matrix (>=1.4-2) class coercion methods

# excursions 2.5.2

* Add `return.marginals.predictor = TRUE` option to tests and examples, for
  new INLA compatibility
* Move inclusion of `omp.h` to before R-related header includes, for clang
  version 13 compatibility

# excursions 2.5.1

* Update links to INLA in documentation

# excursions 2.5.0

* Remove dependency on INLA for tests on CRAN 
* Fix bug for diagonal matrices of size 1x1 
* Add fields URL and BugReports to DESCRIPTION


# excursions 2.4.5

* Fix bug to handle empty sets in continuous interpretations
* Fix compiler warnings; unused variable and mismatching Qinv declaration
* Update COPYRIGHT information

# excursions 2.4.4

* Set fixed buffer size for RngStream.c stream name
* Update Makevars to avoid OpenMP warning

# excursions 2.4.3

* Improved INLA backwards compatibility

# excursions 2.4.2

* Add support for general manifolds in continuous interpretation methods
* Bug fixes to continuous interpretation methods
* Update CITATION information with new JSS manuscript

# excursions 2.4.1

* Minor fixed to C code to avoid warnings.

# excursions 2.4.0

* Added a `NEWS.md` file to track changes to the package.

* Updated repository for INLA

* Add n.iter and seed options for `contourmap()` calculations.

* Add support for the QC method for `contourmap.inla()`.

* Updated documentation.

# excursions 2.3.6

* Fix bug affecting contourmap calculations, caused by the previous fix to `gaussint()`.

# excursions 2.3.5

* Fixed reordering bug in `excursions.variances()` caused by the previous fix to `gaussint()`.

# excursions 2.3.4

* Fixed reordering bug in `gaussint()` that affected standalone use only.

# excursions 2.3.3

* Removed use of deprecated `rBind` and `cBind` and add dependency on `R >= 3.2.0`

# excursions 2.3.2

* Added citation information to DESCRIPTION

# excursions 2.3.1

* Add whitespace after `-f` for make, for wider compatibility, see
  http://pubs/opegroup.org/onlinepubs/9699919799/utilities/make.html

* Updated CITATION information

# excursions 2.3.0
