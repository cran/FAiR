\name{S3methodsFAiR}
\alias{S3methodsFAiR}
\alias{deviance,FA-method}
\alias{df.residual,restrictions-method}
\alias{df.residual,FA-method}
\alias{fitted,FA-method}
\alias{fitted,restrictions-method}
\alias{influence,FA-method}
\alias{model.matrix,FA-method}
\alias{pairs,FA-method}
\alias{residuals,FA-method}
\alias{rstandard,FA-method}
\alias{simulate,FA-method}
\alias{weights,FA-method}
\title{S3 methods for "FA" objects}
\description{
  These methods are technically S4 but are the result of making well-known S3 methods
  into S4 generic functions and defining methods for objects of \code{\link{FA-class}}
  and / or \code{\link{restrictions-class}}. In any event, they provide somewhat standard
  post-estimation functions for factor analysis models.
}
\usage{
\S4method{deviance}{FA}(object)
\S4method{df.residual}{FA}(object)
\S4method{df.residual}{restrictions}(object)
\S4method{fitted}{FA}(object, reduced = TRUE, standardized = TRUE)
\S4method{fitted}{restrictions}(object)
\S4method{influence}{FA}(model)
\S4method{model.matrix}{FA}(object, standardized = TRUE)
\S4method{pairs}{FA}(x, ...)
\S4method{residuals}{FA}(object, standardized = TRUE)
\S4method{rstandard}{FA}(model)
\S4method{simulate}{FA}(object, nsim = 1, seed = NULL, standardized = TRUE)
\S4method{weights}{FA}(object)
}
\arguments{
  \item{object}{An object of \code{\link{FA-class}} or \code{\link{restrictions-class}},
     as appropriate}
  \item{model}{An object of \code{\link{FA-class}}}
  \item{x}{An object of \code{\link{FA-class}}}
  \item{reduced}{Logical indicating whether communalities should be on the
    diagonal of the fitted matrix.}
  \item{standardized}{Logical indicating whether the matrix should be
    calibrated for standardized variables.}
  \item{nsim}{number of simulations}
  \item{seed}{seed to use for random number generation; if \code{NULL} the current
    seed is used}
  \item{\dots}{additional argument(s) for methods}
}
\details{
  The code for each of these methods is quite short. There are some other method 
  definitions for objects that inherit from \code{\link{restrictions-class}} and
  \code{\link{FA-class}} but they differ only in implementation and not in their
  nature or their options.
}
\value{
  \item{deviance}{returns the value of the discrepancy function}
  \item{df.residual}{returns the degrees of freedom}
  \item{fitted}{returns the model's estimate of the covariance or
    correlation matrix among manifest variables in common factor space}
  \item{influence}{returns a square matrix that is equal to 
    \code{residuals() * weights()}}
  \item{model.matrix}{returns the sample covariance or correlation matrix 
    among outcomes}
  \item{pairs}{returns nothing but plots the estimated reference  
    structure correlations in a nice form}
  \item{residuals}{returns a square matrix that contains the 
    difference between \code{model.matrix()} and \code{fitted()} and 
    has uniquenesses along the diagonal (based on correlations by default)}
  \item{rstandard}{covariance residuals standardized by the standard deviations
     of the manifest variables}
  \item{weights}{returns a square matrix with the weights used in the
    discrepancy function. For Yates' weighted least squares estimator
    these weights are as defined in equation 188. For maximum likelihood
    estimation, these weights are proportional to the reciprocal of the 
    crossproduct of the uniquenesses and are only approximately equal to 
    the implied weights that would be used if minimizing the weighted
    sum of squared residuals. For ease of interpretation they are rescaled
    so that the mean weight is 1.0.}
}
\references{
  Yates, A. (1987) \emph{Multivariate Exploratory Data Analysis:
  A Perspective on Exploratory Factor Analysis.} 
  State University of New York Press.
}
\author{Ben Goodrich}
\seealso{\code{\link{loadings}}, \code{\link{cormat}}, and \code{\link{uniquenesses}}}
\examples{
  ## See the example for Factanal()
}
\keyword{ methods }
