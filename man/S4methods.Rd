\name{S4GenericsFAiR}
\docType{genericFunction}
\docType{methods}
\alias{S4GenericsFAiR}
\alias{fitS4}
\alias{bfgs_fitS4}
\alias{gr_fitS4}
\alias{bfgs_helpS4}
\alias{create_start}
\alias{fitS4-methods}
\alias{bfgs_fitS4-methods}
\alias{gr_fitS4-methods}
\alias{bfgs_helpS4-methods}
\alias{create_start-methods}
\alias{fitS4,restrictions,manifest.basic-method}
\alias{bfgs_fitS4,restrictions,manifest.basic-method}
\alias{gr_fitS4,restrictions,manifest.basic-method}
\alias{bfgs_helpS4,restrictions,manifest.basic-method}
\alias{create_start,restrictions,manifest.basic-method}
\title{S4 Generics \& Methods for Package `FAiR'}
\description{
  It is not necessary to understand this help page if one merely
  wants to estimate a factor analysis model. This help page is 
  intended for those who want to modify or extend FAiR or otherwise
  want some idea of how \pkg{FAiR} works \dQuote{behind the scenes}.

  \code{\link{Factanal}} is just a wrapper around a call to 
  \code{\link[rgenoud]{genoud}} in the \pkg{rgenoud} package. These S4
  generic functions are passed to various arguments of
  \code{\link[rgenoud]{genoud}} as part of the optimization process. With
  the exception of the methods for \code{create_start}, it would be very
  unusual for any of these methods to be called directly by the user.
}
\usage{
\S4method{fitS4}{restrictions,manifest.basic}(par, restrictions, manifest, lower, mapping_rule)
\S4method{bfgs_fitS4}{restrictions,manifest.basic}(par, restrictions, manifest, helper, lower)
\S4method{gr_fitS4}{restrictions,manifest.basic}(par, restrictions, manifest, helper, lower)
\S4method{bfgs_helpS4}{restrictions,manifest.basic}(initial, restrictions, manifest, done, lower)
\S4method{create_start}{restrictions,manifest.basic}(number, start, restrictions, manifest, reject)
}
\arguments{
  \item{par}{A numeric vector containing values for all the free parameters 
    to be estimated, which corresponds to the \code{par} argument for 
    \code{\link[rgenoud]{genoud}} and for \code{\link[stats]{optim}}.}
  \item{initial}{Same as \code{par}.}
  \item{restrictions}{An object of \code{\link{restrictions-class}}.}
  \item{manifest}{An object of \code{\link{manifest-class}}.}
  \item{lower}{A small numeric scalar indicating the lower bound for positive
     definiteness or minimum uniqueness; see the corresponding argument
     to \code{\link{Factanal}}.}
  \item{mapping\_rule}{a logical indicating whether to invoke the mapping rule in 
     semi-exploratory models. It takes the value of \code{FALSE} when called by
     \code{\link{optim}}.}
  \item{helper}{A list that contains necessary information for the \code{bfgs_fitS4}
    and \code{gr_fitS4} methods that is passed to the \code{BFGShelp} argument for 
    \code{\link[rgenoud]{genoud}}.}
  \item{done}{A logical indicating whether the optimization has terminated
    or could terminate immediately if a stopping condition is met; see the
    documentation of the \code{BFGShelp} argument to \code{\link[rgenoud]{genoud}}.}
  \item{number}{An integer representing the number of starting values to create}
  \item{start}{A numeric vector or matrix containing initial communality estimates}
  \item{reject}{Logical indicating whether to reject starting values that fail the 
    constraints required by the model}
  \item{\dots}{Further arguments to be passed to downstream functions; not
    currently used.}
}
\details{
  The \code{fitS4} method is responsible for producing a complete numeric \emph{vector} of 
  fit criteria and is passed to the \code{fn} argument of \code{\link[rgenoud]{genoud}}. 
  The default method simply calls \code{\link{restrictions2model}} and then calls an internal
  function that evaluates all the functions in the \code{criteria} slot of the object of 
  \code{\link{restrictions-class}}.

  The \code{bfgs_fitS4} method produces a \emph{scalar} fit criterion, and this 
  method is passed to the \code{BFGSfn} argument of \code{\link[rgenoud]{genoud}}
  and is in turned passed to the \code{fn} argument of \code{\link{optim}}. Usually,
  this criterion is the last criterion produced by the \code{fitS4} method, which is the
  discrepancy function. The \code{gr_fitS4} method produces the gradient
  of the function defined by the \code{bfgs_fitS4} method under traditional maximum
  likelihood estimation; otherwise it is an internal function. Either way, this function
  is passed to the \code{gr} argument of \code{\link[rgenoud]{genoud}} and 
  \code{\link{optim}}. Both \code{bfgs_fitS4} and \code{gr_fitS4} take an
  argument called \code{helper}, which is produced by the \code{bfgs_helpS4}
  method and corresponds to the \code{BFGShelp} argument of 
  \code{\link[rgenoud]{genoud}}. The \code{BFGShelp} method produces a list
  that contains information about the initial value of the genetic individual
  when \code{optim} is called. The \code{bfgs_fitS4} methods and \code{gr_fitS4}
  methods can behave differently depending on the contents of \code{helper}.\cr

  The \code{create_start} method creates a matrix of starting values that is
  then passed to the \code{starting.values} argument of 
  \code{\link[rgenoud]{genoud}}. If you think the starting values are inadequate
  in a particular situation, it is much easier to create a matrix of starting
  values in the global environment and pass it through the \code{\dots} of
  \code{\link{Factanal}} to the \code{starting.values} argument of 
  \code{\link[rgenoud]{genoud}} yourself. See \code{\link{Factanal}} for details
  on doing so.
}
\section{Methods}{
  There are methods for every flavor of \code{\link{restrictions-class}}, although in 
  many cases the methods are simply inherited. Currently, there are only methods for
  objects of class \code{"manifest.basic"}, which are inherited by objects of class
  \code{"manifest.data"} and \code{"manifest.data.mcd"}. However, it is possible to
  tailor methods for different flavors of \code{\link{manifest-class}} or to create 
  new classes that inherit from  \code{\link{restrictions-class}} and write (some) 
  methods for them.
}
\value{
  \item{fitS4}{produces a numeric vector of fits}
  \item{bfgs_fitS4}{produces a numeric scalar as a fit}
  \item{gr_fitS4}{produces a numeric vector that is the gradient at the scalar fit}
  \item{bfgs_helpS4}{produces a list of ancillary material for \code{bfgs_fitS4} and
     \code{gr_fitS4}}
  \item{create_start}{produces a numeric matrix that constitutes the starting population}
}
\author{Ben Goodrich}
\seealso{\code{\link{Factanal}}, \code{\link{manifest-class}}, \code{\link{restrictions-class}} and 
  \code{\link[rgenoud]{genoud}}}
\examples{
  showMethods("fitS4")
  showMethods("bfgs_fitS4")
  showMethods("gr_fitS4")
  showMethods("bfgs_helpS4")
  showMethods("create_start")
}
\keyword{ methods }
\keyword{ multivariate }
\keyword{ models }
