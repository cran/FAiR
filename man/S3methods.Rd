\name{S3methodsFAiR}
\alias{S3methodsFA}
\alias{deviance.FA}
\alias{deviance.FA.general}
\alias{deviance.FA.2ndorder}
\alias{df.residual.FA}
\alias{df.residual.FA.general}
\alias{df.residual.FA.2ndorder}
\alias{fitted.FA}
\alias{fitted.FA.general}
\alias{fitted.FA.2ndorder}
\alias{influence.FA}
\alias{influence.FA.general}
\alias{influence.FA.2ndorder}
\alias{model.matrix.FA}
\alias{model.matrix.FA.general}
\alias{model.matrix.FA.2ndorder}
\alias{pairs.FA}
\alias{pairs.FA.general}
\alias{pairs.FA.2ndorder}
\alias{predict.FA}
\alias{predict.FA.general}
\alias{predict.FA.2ndorder}
\alias{residuals.FA}
\alias{residuals.FA.general}
\alias{residuals.FA.2ndorder}
\alias{rstandard.FA}
\alias{rstandard.FA.general}
\alias{rstandard.FA.2ndorder}
\alias{weights.FA}
\alias{weights.FA.general}
\alias{weights.FA.2ndorder}

\title{S3 methods for objects of class "FA" }
\description{
  These S3 methods for objects of class "FA" or that inherit
  from class "FA" provide fairly standard post-estimation
  functions for factor analysis models.
}
\usage{
\method{deviance}{FA}(object, ...)
\method{deviance}{FA.general}(object, ...)
\method{deviance}{FA.2ndorder}(object, ...)
\method{df.residual}{FA}(object, ...)
\method{df.residual}{FA.general}(object, ...)
\method{df.residual}{FA.2ndorder}(object, ...)
\method{fitted}{FA}(object, ...)
\method{fitted}{FA.general}(object, ...)
\method{fitted}{FA.2ndorder}(object, ...)
\method{influence}{FA}(model, ...)
\method{influence}{FA.general}(model, ...)
\method{influence}{FA.2ndorder}(model, ...)
\method{model.matrix}{FA}(object, ...)
\method{model.matrix}{FA.general}(object, ...)
\method{model.matrix}{FA.2ndorder}(object, ...)
\method{pairs}{FA}(x, ...)
\method{pairs}{FA.general}(x, ...)
\method{pairs}{FA.2ndorder}(x, ...)
\method{predict}{FA}(object, ...)
\method{predict}{FA.general}(object, ...)
\method{predict}{FA.2ndorder}(object, ...)
\method{residuals}{FA}(object, ...)
\method{residuals}{FA.general}(object, ...)
\method{residuals}{FA.2ndorder}(object, ...)
\method{rstandard}{FA}(model, ...)
\method{rstandard}{FA.general}(model, ...)
\method{rstandard}{FA.2ndorder}(model, ...)
\method{weights}{FA}(object, ...)
\method{weights}{FA.general}(object, ...)
\method{weights}{FA.2ndorder}(object, ...)
}
\arguments{
  \item{object}{An object of class 'FA' or that inherits from 'FA'}
  \item{model}{An object of class 'FA' or that inherits from 'FA'}
  \item{x}{An object of class 'FA' or that inherits from 'FA'}
  \item{\dots}{additional argument(s) for methods}
}
\details{
  The code for all of these methods is quite short. Feel free to
  simply call \code{body()} on the \code{method.FA} to see what
  the function does. There is no difference in functionality for
  the methods that inherit from class 'FA' relative to those that
  are defined for class 'FA'.
}
\value{
  \item{deviance}{returns a scalar indicating the discrepancy from
    a perfect fitting model.}
  \item{df.residual}{returns an integer indicating the difference
    between the number of nonredundant elements in the sample correlation
    matrix and the number of estimated parameters.}
  \item{fitted}{returns a square matrix of estimated correlations
    in the common factor space with communalities along the diagonal.}
  \item{influence}{returns a square matrix that is equal to 
    \code{residuals() * weights()}.}
  \item{model.matrix}{returns the sample correlation matrix among
    outcomes.}
  \item{pairs}{returns nothing but plots the estimated reference  
    structure correlations; to repeat: these are the reference structure  
    correlations, rather than the primary pattern coefficients.}
  \item{predict}{returns a matrix with the predicted values of the
    outcomes, which is defined as the product of the factor score
    matrix and the primary pattern coefficient matrix. Thus, an error
    is returned if factor scores were not calculated by 
    \code{\link{Factanal}}.}
  \item{residuals}{returns a square matrix that contains the 
    difference between \code{model.matrix()} and \code{fitted()} and thus 
    has uniquenesses along the diagonal.}
  \item{rstandard}{returns \code{residuals} rescaled into a correlation
    matrix and thus has ones along the diagonal.}
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
\author{Ben Goodrich \url{http://wiki.r-project.org/rwiki/doku.php?id=packages:cran:fair}}
\seealso{\code{\link{confint}}, \code{\link{deviance}}, 
  \code{\link{df.residual}}, \code{\link{fitted}}, \code{\link{influence}},
  \code{\link{model.matrix}}, \code{\link{pairs}}, \code{\link{predict}},
  \code{\link{residuals}}, \code{\link{rstandard}}, \link{FA-class}.}
\examples{
  ## See the example for Factanal()
}
\keyword{ methods }
