# excursions (Development version)

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
