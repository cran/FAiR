\name{FAiR-package}
\alias{FAiR-package}
\alias{FAiR}
\docType{package}
\title{
Factor Analysis in R
}
\description{
This package estimates exploratory, confirmatory, and semi-exploratory factor
analysis models. It utilizes a genetic algorithm, \code{\link[rgenoud]{genoud}} 
for optimization, which makes semi-exploratory factor analysis possible and 
permits new (and enhances old) techniques for transforming factors that are 
extracted via exploratory factor analysis.
}
\details{
\tabular{ll}{
Package: \tab FAiR\cr
Type: \tab Package\cr
Version: \tab 0.0-0\cr
Date: \tab 2008-01-18\cr
License: \tab General Public License version 3 or later\cr
}
The two most important functions are \code{\link{Factanal}} and, in the case
of exploratory factor analysis, \code{\link{Rotate}}. \code{\link{Factanal}}
can also be used to estimate semi-exploratory and confirmatory factor analysis
models, in which case it is unnecessary and impossible to transform the factors.
There are a variety of S4 and S3 methods to facilitate the interpretation of
the results of \code{\link{Factanal}} and \code{\link{Rotate}}. The vignette
has additional information regarding the pop-up menus; execute 
\code{vignette("FAiR")} to read it.
}
\author{
Author and Maintainer: Ben Goodrich
\url{http://wiki.r-project.org/rwiki/doku.php?id=packages:cran:fair}
}

\keyword{ package }
\keyword{ multivariate }
\keyword{ models }
\examples{
  ## See the example for Factanal()
}
