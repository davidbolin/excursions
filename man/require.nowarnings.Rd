\name{require.nowarnings}
\alias{require.nowarnings}
\title{
  Warnings free loading of add-on packages
}
\description{Load and attach add-on packages without warning about
  unavailable packages.
}
\usage{
require.nowarnings(package, lib.loc = NULL, character.only = FALSE)
}
\arguments{
  \item{package}{the name of a package, given as a character string.
}
\item{lib.loc}{
  a character vector describing the location of R library trees
  to search through, or \code{NULL}.  The default value of \code{NULL}
  corresponds to all libraries currently known to
  \code{.libPaths()}.  Non-existent library trees are silently
  ignored.
}
\item{character.only}{
  a logical indicating whether \code{package} can be
  assumed to be a character string.
}
}
\details{
  \code{require.nowarnings(package)} acts the same as
  \code{require(package, quietly = TRUE)}, except that all warnings are
  turned off.  In particular, no warning is given if the package is
  unavailable.
}
\value{
  \code{require.nowarnings} returns (invisibly) a logical indicating whether the
     required package is available.
}

\seealso{
\code{\link{require}}
}
\examples{
## This should produce no output:
if (require.nowarnings(nonexistent)) {
  message("Package loaded successfully")
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
