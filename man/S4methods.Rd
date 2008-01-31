\name{S4GenericsFAiR}
\docType{genericFunction}
\docType{methods}
\alias{S4GenericsFAiR}
\alias{fitS4}
\alias{bfgs_fitS4}
\alias{gr_fitS4}
\alias{bfgs_helpS4}
\alias{create_start}
\alias{create_FAobject}
\alias{fitS4-methods}
\alias{fitS4,numeric,restrictions.factanal-method}
\alias{fitS4,numeric,restrictions.orthonormal-method}
\alias{fitS4,numeric,restrictions.1storder-method}
\alias{fitS4,numeric,restrictions.general-method}
\alias{fitS4,numeric,restrictions.2ndorder-method}
\alias{bfgs_fitS4-methods}
\alias{bfgs_fitS4,numeric,restrictions.factanal,ANY-method}
\alias{bfgs_fitS4,numeric,restrictions.orthonormal,ANY-method}
\alias{bfgs_fitS4,numeric,restrictions.1storder,list-method}
\alias{bfgs_fitS4,numeric,restrictions.general,list-method}
\alias{bfgs_fitS4,numeric,restrictions.2ndorder,list-method}
\alias{gr_fitS4-methods}
\alias{gr_fitS4,numeric,restrictions.factanal,ANY-method}
\alias{gr_fitS4,numeric,restrictions.orthonormal,ANY-method}
\alias{gr_fitS4,numeric,restrictions.1storder,list-method}
\alias{gr_fitS4,numeric,restrictions.general,list-method}
\alias{gr_fitS4,numeric,restrictions.2ndorder,list-method}
\alias{bfgs_helpS4-methods}
\alias{bfgs_helpS4,numeric,restrictions,logical-method}
\alias{bfgs_helpS4,numeric,restrictions.orthonormal,logical-method}
\alias{bfgs_helpS4,numeric,restrictions.1storder,logical-method}
\alias{bfgs_helpS4,numeric,restrictions.general,logical-method}
\alias{bfgs_helpS4,numeric,restrictions.2ndorder,logical-method}
\alias{create_start-methods}
\alias{create_start,restrictions-method}
\alias{create_start,restrictions.factanal-method}
\alias{create_start,restrictions.orthonormal-method}
\alias{create_start,restrictions.1storder-method}
\alias{create_start,restrictions.general-method}
\alias{create_start,restrictions.2ndorder-method}
\alias{create_FAobject-methods}
\alias{create_FAobject,restrictions-method}
\alias{create_FAobject,restrictions.factanal-method}
\alias{create_FAobject,restrictions.orthonormal-method}
\alias{create_FAobject,restrictions.1storder-method}
\alias{create_FAobject,restrictions.general-method}
\alias{create_FAobject,restrictions.2ndorder-method}

\title{S4 Generic Functions and Methods for Package `FAiR'}
\description{
  It is not necessary to understand this help page if one merely
  wants to estimate a factor analysis model. This help page is 
  intended for those who want to modify or extend FAiR or otherwise
  want some idea of how FAiR works ``behind the scenes''.\cr
  These generic functions are used to create starting values for,
  estimate, and encapsulate the results of a factor analysis model
  that is estimated by \code{\link{Factanal}}. The various methods
  for these models pertain to  different types of factor analysis models.
}
\usage{
fitS4(par, object, ...)
bfgs_fitS4(par, object, helper, ...)
gr_fitS4(par, object, helper, ...)
bfgs_helpS4(initial, object, done, ...)
create_start(restrictions, ...)
create_FAobject(restrictions, ...)
}
\arguments{
  \item{par}{A numeric vector containing values for all the free or 
    parameters to be estimated, which corresponds to the \code{par} argument 
    for \code{\link[rgenoud]{genoud}} and for \code{\link[stats]{optim}}.}
  \item{initial}{Same as \code{par}.}
  \item{object}{In FAiR, this object is always inherits from class
    \code{"restrictions"}, see \code{\link{restrictions-class}}.}
  \item{helper}{An object, typically a list that contains necessary information 
    for the \code{bfgs_fitS4} and \code{gr_fitS4} methods as documented in
    \code{\link[rgenoud]{genoud}}.}
  \item{done}{A logical indicating whether the optimization has terminated
    or could terminate immediately if a stopping condition is met.}
  \item{restrictions}{An object of that inherits from class \code{"restrictions"},
    see \code{\link{restrictions-class}}.}
  \item{\dots}{Further arguments to be passed to downstream functions}
}
\details{
  None of these S4 generic functions have default methods; any additional
  methods must be tailored to the specific factor analysis model being 
  estimated.

  The \code{fitS4} method is the function for which the optimal \code{par} is
  sought and produces a numeric \emph{vector} of fit criteria. This method is passed 
  to the \code{fn} argument of \code{\link[rgenoud]{genoud}} for \emph{lexical}
  optimization. Thus, the last element of this vector should be fully continuous 
  so that all ties are broken among unique individuals in the population of the 
  genetic algorithm. For example, the last element is the value of the 
  log-likelihood function when \code{method = "MLE"} is specified in the call to 
  \code{\link{Factanal}}.

  Lexical optimization can be seen as optimization with respect to the last
  criterion produced by the \code{fitS4} method among indiivduals in the population
  of the genetic algorithm with the same values on all previous criteria produced
  by the \code{fitS4} method. Thus, previous elements of this vector of criteria
  are values of \dQuote{constraint} criteria and by convention are operationalized
  as piecewise functions of \code{par} that take the value of 1.0 if the constraint 
  is satisfied and some value less than 1.0 when the constraint is not satisfied.
  For example, the constraint that prevents Heywood cases takes the value of 1.0
  when all specific variances are positive and takes the value of the minimum 
  specific variance otherwise. In FAiR, \code{fitS4} is supplied with an additional
  argument, \code{S}, which is the sample correlation matrix among outcomes.

  The \code{bfgs_fitS4} method produces a \emph{scalar} fit criterion, and this 
  method is passed to the \code{BFGSfn} argument of \code{\link[rgenoud]{genoud}}
  and is in turned passed to the \code{fn} argument of \code{\link{optim}}. Usually,
  this criterion is the last criterion produced by the \code{fitS4} method, 
  such as the log-likelihood. The \code{gr_fitS4} method produces the gradient
  of the function defined by the \code{bfgs_fitS4} method, and this method is
  passed to the \code{gr} argument of \code{\link[rgenoud]{genoud}} and 
  \code{\link{optim}}. Both \code{bfgs_fitS4} and \code{gr_fitS4} take an
  argument called \code{helper}, which is produced by the \code{bfgs_helpS4}
  method and corresponds to the \code{BFGShelp} argument of 
  \code{\link[rgenoud]{genoud}}. The \code{BFGShelp} method typically produces
  a list that contains information about the initial value of the individual
  when \code{optim} is called. In FAiR, the \code{bfgs_fitS4} methods and
  \code{gr_fitS4} methods behave differently depending on the contents of
  \code{helper}. In FAiR, each of these methods are supplied with an additional
  argument, \code{S}, which is the sample correlation matrix among outcomes.\cr

  The \code{create_start} method creates a matrix of starting values that is
  then passed to the \code{starting.values} argument of 
  \code{\link[rgenoud]{genoud}}. If you think the starting values are inadequate
  in a particular situation, it is much easier to create a matrix of starting
  values in the global environment and pass it through the \code{\dots} of
  \code{\link{Factanal}} to the \code{starting.values} argument of 
  \code{\link[rgenoud]{genoud}} yourself. See \code{\link{Factanal}} for details
  on doing so. It is only necessary to write a \code{create_start} method
  when extending the \code{\link{restrictions-class}} to implement a new kind
  of factor analysis model. In FAiR, the \code{create_start} method is
  supplied with additional arguments: \code{pop.size}, which corresponds to
  the argument of the same name in \code{\link[rgenoud]{genoud}} and indicates
  how many starting vectors to create; \code{start} which is a one-column matrix
  containing the initially estimated communalities for the outcomes; and \code{S},
  which is the sample correlation matrix among the outcome variables.

  Finally, the \code{create_FAobject} takes the output of 
  \code{\link[rgenoud]{genoud}} and produces an object of class "FA" or an
  object that inherits from class "FA". In FAiR, the \code{create_FAobject}
  method is supplied with additional arguments: \code{opt}, which is the 
  list produced by \code{\link[rgenoud]{genoud}}; \code{manifest}, which
  is a list containing the sample covariance matrix, sample correlation matrix,
  and the number of observations in the sampe; and \code{call}, which is the
  call to \code{\link{Factanal}}.

  To summarize, in order to write a new method for FAiR, create a new class that
  inherits from \code{\link{restrictions-class}}. Write a \code{fitS4} method
  for this new class that calculates all criteria that you want to find the
  \emph{lexical} optimum of, given \code{par}, this object, and \code{S}. Write a 
  \code{bfgs_fitS4} method to calculate some scalar criterion. The \code{bfgs_fitS4} 
  method would typically reuse a lot of code from the \code{fitS4} method or, in the 
  case of all \code{bfgs_fitS4} methods currently implemented in FAiR, simply calls 
  the \code{fit_S4} method and selects the first value of the vector that is not
  equal to \eqn{1.0}. Write a \code{gr_fitS4} method that calculates the gradient
  of the \code{bfgs_fitS4} method at \code{par}, which will also reuse a lot of
  code from the \code{fitS4} method. If you want to do anything fancy with the
  \code{bfgs_fitS4} and \code{gr_fitS4} methods write a \code{bfgs_helpS4} method 
  that evaluates the \code{fitS4} method at \code{initial} and passes the result
  to the \code{bfgs_fitS4} and \code{gr_fitS4} methods via their \code{helper}
  argument. Finally, write a \code{create_FAobject} method that produces an
  object of class \code{"FA"} or that inherits from class \code{"FA"} (in which
  case you would need to define that class and at least write \code{summary} 
  methods for it). It is recommended to also write a \code{create_start} method
  that creates starting values in light of the \code{fitS4} method. Starting
  values should ideally be widely dispersed throughout the subset of the 
  parameter space where all of the constraints are satisfied. For each of these
  methods, the existing methods in FAiR provide a good template to start from.
}
\section{Methods}{
  See the examples below. The essential concept is that the estimation behavior
  depends on \code{object} and \code{restrictions} which is an object of
  a class that inherits from the \code{"restrictions"} superclass. See also 
  \code{\link{restrictions-class}}.
}
\value{
  \item{fitS4}{produces a numeric vector}
  \item{bfgs_fitS4}{produces a numeric scalar}
  \item{gr_fitS4}{produces a numeric vector}
  \item{bfgs_helpS4}{typically produces a list}
  \item{create_start}{produces a numeric matrix}
  \item{create_FAobject}{produces an object that is or inherits from
    class "FA"}
}
\author{Ben Goodrich \url{http://wiki.r-project.org/rwiki/doku.php?id=packages:cran:fair}}
\seealso{\code{\link{Factanal}}, \code{\link{restrictions-class}} and 
  \code{\link[rgenoud]{genoud}}}
\examples{
  showMethods("fitS4")
  showMethods("bfgs_fitS4")
  showMethods("gr_fitS4")
  showMethods("bfgs_helpS4")
  showMethods("create_start")
  showMethods("create_FAobject")
}
\keyword{ methods }
\keyword{ multivariate }
\keyword{ models }
