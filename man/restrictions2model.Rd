\name{restrictions2model}
\docType{methods}
\alias{restrictions2model}
\alias{restrictions2model-methods}
\alias{restrictions2model,restrictions,manifest-method}
\title{Map free parameters to matrices}
\description{
  It is not necessary to understand this help page if one merely
  wants to estimate a factor analysis model. This help page is 
  intended for those who want to modify or extend FAiR or otherwise
  want some idea of how FAiR works \dQuote{behind the scenes}.

  This S4 generic function takes a numeric vector of free parameters and
  manipulates it into the matrices that are typical when estimating a
  factor analysis model. There is little need to call it directly. Various 
  methods are defined for the \code{restrictions} object, corresponding to 
  different kinds of factor analysis models; see \code{\link{FAiR-package}}.
}
\section{Methods}{
  The arguments that are part of the signature are \code{restrictions} and
  \code{manifest}. Methods are currently only defined for objects of class
  \code{"manifest.basic"}, which are inherited by objects of class 
  \code{"manifest.data"} and \code{"manifest.data.mcd"}. There are methods
  for each of the classes that inherit from \code{\link{restrictions-class}},
  except for \code{"restrictions.factanal"}, which does not utilize the
  \code{restrictions2model} mechanism.

  There are also two arguments that are \strong{not} part of the signature.
  The first is \code{par}, which is a numeric vector of free parameters and
  is conceptually similar to the \code{par} argument to \code{\link{optim}}.
  The second is \code{lower}, which is a small positive number that is used
  as a threshold for positive-definiteness of various matrices.
}
\usage{
\S4method{restrictions2model}{restrictions,manifest}(par, restrictions, manifest, lower, mapping_rule)
}
\arguments{
  \item{par}{A numeric vector containing values for all the free parameters 
    to be estimated, which corresponds to the \code{par} argument for 
    \code{\link[rgenoud]{genoud}} and for \code{\link[stats]{optim}}.}
  \item{restrictions}{An object of \code{\link{restrictions-class}}.}
  \item{manifest}{An object of \code{\link{manifest-class}}.}
  \item{lower}{A small numeric scalar indicating the lower bound for positive
     definiteness or minimum uniqueness; see the corresponding argument
     to \code{\link{Factanal}}.}
  \item{mapping\_rule}{a logical that controls whether the mapping rule is invoked in
     semi-exploratory factor analysis. This argument is used to circumvent
     the mapping rule in some circumstances, like when the function is called
     indirectly by \code{\link{optim}}}
}
\details{
This function is called internally by \code{\link{Factanal}} thousands of times
during the course of the optimization. Let the factor analysis model in the population be
\deqn{\Sigma = \Omega(\beta\Phi\beta^\prime + \Theta)\Omega}{Sigma = Omega(beta Phi beta' + Theta)Omega}
and let \eqn{\theta}{theta} be a vector of all the free parameters in the factor
analysis model. The \code{restrictions2model} methods are essentially a mapping 
from \eqn{\theta}{theta} to the free elements of \eqn{\beta}{beta}, 
\eqn{\Phi}{Phi}, and \eqn{\Omega}{Omega}. The methods are currently defined on
the \code{restrictions} argument, which is fairly skeletal at the outset, is
filled in by the body of the \code{restrictions2model} method using a bunch of
calls to the \code{\link{make_parameter-methods}}.

The \code{restrictions2model} method can (and does) also do some checking of
whether the parameters in \code{par} are admissable in the context of a factor
analysis model. For each admissability check, if \code{par} is admissable, it
should receive the value of \eqn{-1.0} in a numeric vector whose length is equal
to the number of admissability checks. If \code{par} is inadmissable with respect
to some check, it should receive a value greater than \eqn{-1.0}. This numeric
vector should also be returned as an element of a list called \code{"criteria"}.
}
\seealso{\code{\link{make_parameter}}, \code{\link{restrictions-class}}}
\author{Ben Goodrich}
\value{
A list with exactly two elements named
\item{criteria}{a numeric vector}
\item{restrictions}{a filled-in version of the \code{restrictions} object}
}
\examples{
  showMethods("restrictions2model")
}

\keyword{methods}
\keyword{multivariate}
