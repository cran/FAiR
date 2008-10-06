\name{loadings}
\docType{genericFunction}
\docType{methods}
\alias{loadings}
\alias{cormat}
\alias{uniquenesses}
\alias{loadings-methods}
\alias{cormat-methods}
\alias{uniquenesses-methods}
\alias{coef,FA-method}
\alias{coef,restrictions-method}
\alias{coef,parameter.coef-method}
\alias{cormat,FA-method}
\alias{cormat,FA.2ndorder-method}
\alias{cormat,restrictions-method}
\alias{cormat,restrictions.2ndorder-method}
\alias{cormat,parameter.cormat-method}
\alias{loadings,ANY-method}
\alias{loadings,FA-method}
\alias{loadings,FA.general-method}
\alias{loadings,restrictions-method}
\alias{loadings,restrictions.general-method}
\alias{uniquenesses,FA-method}
\alias{uniquenesses,FA.general-method}
\alias{uniquenesses,restrictions-method}
\title{Extractor ("get") functions}
\description{
These functions extract various elements of formal S4 objects that are
important in factor analysis models, namely the loadings, the correlations
among factors, and the unique variances. Occasionally, it may be useful to
call these generic functions directly.
}
\section{Methods}{
There are methods for every flavor of \code{\link{FA-class}}
and virtually all flavors of \code{\link{restrictions-class}}. Also,
in the code of \code{cormat}, there is a method for objects that inherit
from \code{\link{manifest-class}}.
}
\usage{
\S4method{coef}{FA}(object)
\S4method{coef}{restrictions}(object)
\S4method{loadings}{FA}(x, matrix = "PP", standardized = TRUE)
\S4method{loadings}{FA.general}(x, matrix = "PP", standardized = TRUE, level = 1)
\S4method{loadings}{restrictions.general}(x, standardized = TRUE, level = 1)
\S4method{cormat}{FA}(object, matrix = "PF")
\S4method{cormat}{FA.2ndorder}(object, matrix = "PF", level = 1)
\S4method{cormat}{restrictions}(object)
\S4method{cormat}{restrictions.2ndorder}(object, level = 1)
\S4method{uniquenesses}{FA}(object, standardized = TRUE)
\S4method{uniquenesses}{FA.general}(object, standardized = TRUE, level = 1)
\S4method{uniquenesses}{restrictions}(object, standardized = TRUE)
}
\arguments{
  \item{object}{an object that inherits from \code{\link{FA-class}} or \code{\link{restrictions-class}}}
  \item{x}{an object that inherits from \code{\link{FA-class}} or \code{\link{restrictions-class}}}
  \item{matrix}{a character string with exactly two letters indicating which matrix to
    extract; see the Details section}
  \item{standardized}{a logical indicating whether to standardize the result so that it
    is calibrated for a correlation matrix among manifest variables, rather than their 
    covariance matrix}
  \item{level}{either \eqn{1} or \eqn{2} to indicate from which level of the factor analysis
    model is pertinent when the model has two levels}
}
\details{
Let the factor analysis model be
\deqn{\Sigma = \Omega(\beta\Phi\beta^\prime + \Theta)\Omega}{Sigma = Omega(beta Phi beta' + Theta)Omega}
By default, the \code{loadings} methods extract the estimate of \eqn{\beta}{beta}, 
the \code{cormat} methods extract the estimate of \eqn{\Phi}{Phi}, and the 
\code{uniquenesses} methods extract the diagonal of \eqn{\Theta}{Theta}. In addition,
the \code{coef} methods and the \code{loadings} methods that are defined for objects
\code{restrictions-class} extract the primary pattern matrix (at level 1). 

At the moment there is no special function to get the diagonal of \eqn{\Omega}{Omega}, which
is a diagonal matrix of \emph{estimated} standard deviations of the manifest variables. However,
they can be extracted from the appropriate slot using the \code{@} operator. Also, if
\code{standardized = FALSE} in the call to \code{loadings} or \code{uniquenesses}, then the 
loadings or uniquenesses are scaled by these estimated standard deviations to produce
estimates on the covariance metric.

Additionally, for the \code{loadings} and \code{cormat} methods that are defined on objects of
\code{\link{FA-class}}, the \code{matrix} argument can be specified to extract a different
set of estimated coefficients or correlations. By default, \code{matrix = "PP"} for these 
\code{loadings} methods, indicating that the primary pattern matrix should be extracted. 
Other possible choices are \code{"PS"} to extract the primary structure matrix (defined as 
\eqn{\beta\Phi}{beta Phi}), \code{"RS"} to extract the reference structure matrix (which is 
column-wise proportional to \eqn{\beta}{beta}), \code{"RP"} to extract the reference pattern
matrix (which is column-wise proportional to \eqn{\beta\Phi}{beta Phi}), and \code{"FC"}
to extract the factor contribution matrix (which is defined as
\eqn{\beta * (\beta\Phi)}{beta * (beta Phi)}, where the \code{*} indicates element-by-element
multiplication of two matrices with the same dimensions).

By default, \code{matrix = "PF"} for these \code{cormat} methods, indicating that the
correlation matrix among primary factors should be extracted. Other possible choices
are \code{"RF"} to extract the correlation matrix among reference factors and 
\code{"PR"} to extract the (diagonal) correlation matrix between primary and reference
factors.

In the case of a two-level model, the \code{level} argument can be specified to 
extract such matrices from the second level of the model (including the methods
for the \code{uniquenesses} generic).
}
\author{Ben Goodrich}
\value{
\code{loadings} outputs a matrix of S3 class \code{"loadings"}, which has a special
print method (see \code{\link{print.loadings}}). \code{coef} returns the primary pattern
matrix at level one and is not of class \code{"loadings"}. The \code{cormat} methods 
output a (symmetric) matrix, and the \code{uniquenesses} methods output a non-negative 
numeric vector.
}
\examples{
  ## See the example for Factanal()
}

\keyword{methods}
