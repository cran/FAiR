\name{FAiR-package}
\alias{FAiR-package}
\alias{FAiR}
\concept{exploratory factor analysis}
\concept{EFA}
\concept{confirmatory factor analysis}
\concept{CFA}
\concept{semi-exploratory factor analysis}
\concept{SEFA}
\docType{package}
\title{
Factor Analysis in R
}
\description{
This package estimates exploratory, confirmatory, and semi-exploratory factor
analysis models via a genetic algorithm, namely \code{\link[rgenoud]{genoud}}.
This use of a genetic algorithm is tantamount to restricted optimization with
virtually unlimited possibilities for restrictions. In particular, semi-exploratory
factor analysis, which is new to the literature, minimizes a discrepancy function
subject to a restriction on the \emph{number} of exact zeros in each column of the 
primary pattern matrix but does not require that the \emph{locations} of the zeros 
be specified in advance (as in confirmatory factor analysis). \pkg{FAiR} encourages
the use of \emph{inequality} restrictions on functions of multiple parameters to 
characterize testable hypotheses.
}
\details{
\tabular{ll}{
Package: \tab FAiR\cr
Type: \tab Package\cr
URL: \tab \url{http://wiki.r-project.org/rwiki/doku.php?id=packages:cran:fair}\cr
Version: \tab 0.4-7\cr
Date: \tab 2010-02-24\cr
License: \tab AGPL (>= 3.0) but see the LICENSE file for details
}

Let the factor analysis model in the population be
\deqn{\Sigma = \Omega(\beta\Phi\beta^\prime + \Theta)\Omega}{Sigma = Omega(beta Phi beta' + Theta)Omega}
where \eqn{\Sigma}{Sigma} is the covariance matrix among outcome variables, \eqn{\Omega}{Omega}
is a diagonal matrix of standard deviations of the manifest variables, \eqn{\beta}{beta} is the 
primary pattern matrix (calibrated to standardized variables) with one column per factor, 
\eqn{\Phi}{Phi} is the correlation matrix among the primary factors, and \eqn{\Theta}{Theta} 
is the diagonal matrix of uniquenesses, which is fully determined by \eqn{\beta}{beta}, 
\eqn{\Phi}{Phi}, and the requirement that the matrix within parentheses has ones down its diagonal. 
Hence, \eqn{\beta\Phi\beta^\prime + \Theta}{beta Phi beta' + Theta} is the model's purported correlation 
matrix among outcome variables as a function of the factors.

Each of the matrices on the right-hand side is a parameter to be estimated, and unlike many structural
equation modeling programs, there is no mechanism for \dQuote{translating} the model from a path
diagram or otherwise avoiding the matrix algebra representation of the model. On a technical
programming note, each of these matrices is represented by a \code{\link{parameter-class}} in \pkg{FAiR}, 
which includes a slot for the (proposed) estimate but also includes slots for ancillary information.

The usual steps to estimate and interpret a factor analysis model are as follows:

0. Get your data into \R somehow. It is best to load the raw data in one of the usual
fashions (e.g. \code{\link{read.table}}, \code{\link[foreign]{read.spss}}, etc.) and the 
\code{\link{read.cefa}} function can also be used if your data are saved in the format 
used by CEFA 2.0. If you only have a covariance matrix, then \code{\link{read.triangular}} 
can be used to load it into \R.

1. Call \code{\link{make_manifest}} to construct the left-hand side of
the factor analysis model, namely an S4 object to house the sample estimate of \eqn{\Sigma}{Sigma}
and some other information (e.g. number of observations).

2. Call \code{\link{make_restrictions}} to establish the additional restrictions to
be imposed on the right-hand side of factor analysis model, inclusive of whether the model
is exploratory, semi-exploratory, or confirmatory and what discrepancy function
to use. There is an extensive GUI that pops up when \code{\link{make_restrictions}} is
called to guide you through this step. \pkg{FAiR} differs fundamentally from other factor 
analysis software in that it permits you to impose inequality restrictions on functions of 
\eqn{\beta}{beta} and \eqn{\Phi}{Phi}. Hence, the \code{\link{restrictions-class}} is critical
to the way \pkg{FAiR} is programmed internally and houses S4 objects representing each of the 
matrices to be estimated.

3. Call \code{\link{Factanal}} to estimate the model and thereby produce estimates of
\eqn{\Omega}{Omega}, \eqn{\beta}{beta}, \eqn{\Phi}{Phi}, and \eqn{\Theta}{Theta}.

4. (exploratory factor analysis only) Call \code{\link{Rotate}} to choose
a transformation matrix (\bold{T}) for the factors. There is an extensive GUI that pops up when
\code{\link{Rotate}} is called to guide you through this step. Again, \code{\link{Rotate}} 
differs fundamentally from other approaches to factor rotation in that it permits you to
impose inequality restrictions on functions of parameters when searching for \bold{T}.

5. Call the usual post-estimation methods to interpret the estimates, like 
\code{\link{summary}}, \code{\link{pairs,FA-method}}, etc., and call
\code{\link{model_comparison}} to see the test statistics and fit indices.

The vignette has additional information regarding the pop-up menus produced in step 2 by\cr 
\code{\link{make_restrictions}} and in step 4 by \code{\link{Rotate}}; execute 
\code{vignette("FAiR")} to read it. The primary examples are in \code{\link{Factanal}} and
\code{\link{Rotate}}.
}
\author{
Ben Goodrich
}
\examples{
  ## See the examples for Factanal() and Rotate()
}
\keyword{ package }
\keyword{ multivariate }
\keyword{ models }
