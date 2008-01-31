\name{make_restrictions}
\alias{make_restrictions}
\title{ Make an object of class "restrictions" }
\description{
  This function is a constructor for objects that are of class 
  \code{"restrictions"} or that inherit from \code{"restrictions"}.
  Rarely would a user call this function directly, since it is
  called by \code{\link{Factanal} in the usual case}.
}
\usage{
make_restrictions(factors, model, method, fixed, covmat, criteria = NULL)
}
\arguments{
\item{factors}{ The number (>0) of factors to be fitted, which differs
    from the argument in \code{\link{factanal}} in that 
    \code{factors} can be a numeric vector of length two to indicate
    the number of factors to extract at level one and level two of a 
    two-level semi-exploratory or confirmatory factor analysis model.
    If a single number, \code{make_restrictions} will ask whether to estimate
    a second level if appropriate. }  
  \item{model}{ A character string indicating "SEFA", "EFA", or "CFA" to
     indicate whether a semi-exploratory, an exploratory, or a confirmatory 
     factor analysis model should be estimated. Defaults to "SEFA". }
  \item{method}{ A character string indicating "MLE" or "YWLS" to
     indicate how the \code{model} should be estimated. Defaults to "MLE".
     The "YWLS" option uses Yates' (1987) weighted-least squares criterion
     as opposed to most of the weighted-least squares criteria that are
     usually mentioned in the factor analysis literature.}
  \item{fixed}{ An optional matrix or list of two matrices that specifies
     the values of certain coefficients. Unrestricted coefficients should
     be denoted by \code{NA}. If \code{fixed} is a matrix, it should have
     rows equal to the number of outcome variables and columns equal to
     the number of factors at level one. If \code{fixed} is a list of two
     two matrices, the first element of the list corresponds to the 
     coefficient matrix at level one and the second element corresponds to the 
     coefficient matrix at level two (and should have rows equal to the 
     number of first-order factors and columns equal to the number of
     second-order factors. By definition, providing \code{fixed} is
     inappropriate if \code{model = "EFA"}, but \code{fixed} need not be
     supplied in any case. If appropriate, \code{make_restrictions} will prompt 
     you about fixing coefficients to certain values with its pop-up menus.}
  \item{covmat}{ A covariance matrix, or a covariance list as returned by
    \code{\link{cov.wt}} or similar. It should have row and column names.}
  \item{criteria}{ An \strong{optional} list whose elements should be functions or 
     character strings that name functions to be used as criteria during
     the lexical optimization when \code{model != "EFA"}. It is almost always
     best to leave this argument unspecified to allow \code{make_restrictions} 
     to prompt for these criteria.\cr

     If \code{criteria} is a list that includes character strings
     the strings should be one or more of \code{"no_suppressors_xxx"},
     \code{"dets_xxx"}, or \code{"cohyperplanarity"} where \code{"xxx"}
     is either \code{"1st"} or \code{"2nd"} to indicate whether the
     criterion should be applied to the first or second level of the model.
     Thus, a suffix of \code{"2nd"} is only appropriate if a two-level
     model is estimated. The function implied by \code{method} is automatically
     appended to the end of this list to serve as the ultimate lexical 
     criterion.}
}
\details{
  This function is essential but it is called internally by \code{\link{Factanal}}
  when the \code{restrictions} argument of \code{\link{Factanal}} is not specified.
  Thus, a typical user would never need to call this function directly. It is 
  somewhat convenient when conducting simulations or debugging, in the sense that 
  the appropriate object can be created once and then passed repeatedly to the
  \code{restrictions} argument of \code{\link{Factanal}} to avoid having to 
  repeatedly answer the questions in the pop-up menus.
}
\value{
  Returns an object of class \code{"restrictions"} or that inherits from
  class \code{"restrictions"}. This object would then be passed to the
  \code{restrictions} argument of \code{\link{Factanal}}.
}
\author{Ben Goodrich \url{http://wiki.r-project.org/rwiki/doku.php?id=packages:cran:fair}}
\seealso{ \code{\link{Factanal}} and \code{\link{restrictions-class}} }
\examples{
  \dontrun{
    res <- make_restrictions(covmat = Harman74.cor) # answer pop-up questions
    show(res)
  }
}
\keyword{ multivariate }
