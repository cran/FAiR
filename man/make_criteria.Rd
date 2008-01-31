\name{make_criteria}
\alias{make_criteria}
\title{ Create Lexical Criteria for Factor Transformation }
\description{
  This function establishes the lexical criteria that are used by
  \code{\link{Rotate}} to find the optimal transformation of the factors
  Rarely would a user call this function directly, since it is
  called by \code{\link{Rotate} in the usual case}.
}
\usage{
make_criteria(factors, weights, FAobject)
}
\arguments{
  \item{factors}{ The number of factors }
  \item{weights}{ This argument is only relevant if the varphi criterion
    is used for transformation and should either be NULL or a numeric vector 
    giving the weights to be used If NULL, then a dynamic weighting scheme 
    will be used for the with the varphi criterion.}
  \item{FAobject}{ An object of class \code{"FA"} produced by \code{\link{Factanal}}}
}
\details{
  This function is essential but it is called internally by \code{\link{Rotate}}
  when the \code{criteria} argument of \code{\link{Rotate}} is not specified.
  Thus, a typical user would never need to call this function directly. It is 
  somewhat convenient when conducting simulations or debugging, in the sense that 
  the appropriate object can be created once and then passed repeatedly to the
  \code{criteria} argument of \code{\link{Rotate}} to avoid having to 
  repeatedly answer the questions in the pop-up menus.

}
\value{
  Returns a list of functions of length at least one
}
\author{Ben Goodrich \url{http://wiki.r-project.org/rwiki/doku.php?id=packages:cran:fair}}
\seealso{ \code{\link{Rotate}} }
\examples{
  \dontrun{
    criteria <- make_criteria(factors = 3, weights = NULL) # answer pop-ups
  }
}
\keyword{ multivariate }
