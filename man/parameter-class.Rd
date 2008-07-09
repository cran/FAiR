\name{parameter-class}
\docType{class}
\docType{genericFunction}
\docType{methods}
\alias{parameter-class}
\alias{parameter.cormat-class}
\alias{parameter.scale-class}
\alias{parameter.coef-class}
\alias{parameter.coef.nl-class}
\alias{parameter.coef.SEFA-class}
\alias{parameter.coef.SEFA.nl-class}
\alias{make_parameter}
\alias{make_parameter-methods}
\alias{show,parameter.cormat-method}
\alias{show,parameter.coef-method}
\title{Class "parameter"}
\description{
  It is not necessary to understand this help page if one merely
  wants to estimate a factor analysis model. This help page is 
  intended for those who want to modify or extend FAiR or otherwise
  want some idea of how \pkg{FAiR} works \dQuote{behind the scenes}.

  This class is used internally to hold information about the 
  parameters in a factor analysis model.
}

\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("parameter", ...)}.
However, rarely, if ever, would one want to do so directly because
\code{\link{make_restrictions}} instantiates them and puts them
into the slots of an object of \code{\link{restrictions-class}}.
}
\section{Slots}{
  The \code{"parameter"} class is virtual and has the following slots:
  \describe{
    \item{\code{x}:}{Object of class \code{"ANY"} but is a
      numeric vector or numeric matrix in all known inherited
      classes. Its elements can be free parameters or parameters
      that are fixed and not estimated. In the case of the former,
      the corresponding element of \code{x} should be \code{NA}.
      In the case of the latter, the corresponding element of 
      \code{x} should be the number that the parameter is fixed to
      or \code{Inf} if the parameter is \dQuote{fixed} to be a
      function of other parameters. }
    \item{\code{free}:}{Object of class \code{"ANY"} but for all 
      known inherited class is a logical vector or logical matrix whose 
      length is the same as the length of \code{x} that indicates which
      elements of \code{x} are considered \dQuote{free} parameters.}
    \item{\code{num_free}:}{Integer indicating the number of free
      parameters in \code{x}.}
    \item{\code{select}:}{A logical vector indicating which elements
      of the grand parameter vector correspond to the free elements
      of \code{x}.}
    \item{\code{invalid}:}{A numeric scalar that is \eqn{0.0} if
      the parameters are \emph{not} invalid and some number other
      than \eqn{0.0} if the parameters are invalid.}
    \item{\code{Domains}:}{an array that is either empty or stacks exactly
      two matrices that are the same dimension as \code{x}. The lower
      matrix should indicate the lower bound for the parameter and
      the upper matrix should indicate the upper bound. It is often
      unnecessary to specify \code{Domains} because they are completed
      by the \code{\link{make_restrictions-methods}}.}
  }
  The \code{"parameter.cormat"} class extends \code{"parameter"} and
  is used for correlation matrices among primary factors. Hence,
  \code{x} is of class \code{"matrix"} and the only free parameters
  are in the lower triangle; however, not all elements of the lower
  triangle need to be free parameters. The \code{"parameter.scale"}
  class extends \code{"parameter"} and is used for (diagonal) matrices
  of standard deviations of the outcome variables that pre- and 
  post-multiply the correlation matrix in the \dQuote{embedded correlation}
  parameterized currently used for all models in \pkg{FAiR}. Hence, 
  \code{x} is of class \code{"numeric"} and is the diagonal of such a 
  scale matrix.Typically, all diagonal elements are free but not always, 
  as in the case where some observed variable is designated as a factor,
  in which case its standard deviation is estimated from the data
  rather than by the model.

  The \code{"parameter.coef"} class extends \code{"parameter"} and is used 
  for primary pattern matrices. Hence, \code{"x"} is of class (numeric)
  \code{"matrix"} and its cells may be free (designated by \code{NA}) or
  fixed. It has one additional slot,
  \describe{
    \item{\code{equalities}:}{a list (possibly with zero length) of objects 
      of \code{\link{equality_restriction-class}} to indicate equality restrictions
      among some cells of \code{x}.}
  }
  The \code{"parameter.coef.nl"} class extends \code{"parameter.coef"} and has 
  one additional slot,
  \describe{
    \item{\code{nonlinearities}:}{a function with an argument called \code{"pattern"}. 
      After \code{x} is filled with free parameters and any equality restrictions are 
      resolved, \code{x} is passed to this function whose body can enforce other exact 
      restrictions on the cells. For example, one \dQuote{"fixed"} cell could be the 
      product of other cells. Then, this function must return this primary pattern 
      matrix whose cells are all finite numbers, which will be reinserted into the
      \code{x} slot for later use.}
  }

  The \code{"parameter.coef.SEFA"} class extends \code{"parameter.coef"} and
  is used in semi-exploratory factor analysis. It has two additional slots,
  \describe{
    \item{\code{rankcheck}:}{either \code{"reiersol"} or \code{"howe"} to indicate
      which theorem is to be used in checking the rank of submatrices of the
      primary pattern matrix with exact zeros in a column. See the Details section.}
    \item{\code{mapping_rule}:}{a function that defaults to the one documented in
      \code{\link{mapping_rule}} and is used to make some cells of \code{x}
      exactly zero. It is possible to define a different function for this slot
      but see \code{\link{mapping_rule}} because the underlying code is 
      hard-coded to assume some of the same arguments in some places. This function
      is called after the cells of \code{x} have been filled and any equality
      restrictions have been resolved and should output a matrix that will be 
      reinserted into the \code{x} slot.}
    \item{\code{squashed}:}{a logical matrix indicating which cells of the coefficient
      matrix were squashed to zero by the mapping rule. This slot is typically empty
      until the model has been estimated, at which point the 
      \code{\link{create_FAobject-methods}} must fill it.}
  }
  The \code{"parameter.coef.SEFA.nl"} class extends  \code{"parameter.coef.SEFA"}
  and has one additional slot,
  \describe{
    \item{\code{nonlinearities}:}{A function, the same as above for the
      \code{"parameter.coef.nl"} class. The function enforcing nonlinear 
      restrictions is called before the mapping rule function.}
  }
}
\section{Methods}{
  A \code{\link{cormat}} method extracts the \code{x} slot of an object of
  \code{"parameter.cormat"}. A \code{\link{coef}} method extracts the 
  \code{x} slot of an object that inherits from \code{"parameter.coef"}.
  There is no special extractor function for objects of \code{"parameter.scale"}.
  There are also \code{\link{show}} methods.

  The most important methods are those defined for the \code{make_parameter}
  S4 generic. This function takes two arguments, \code{"par"}, which is a
  numeric vector but not part of the signature, and \code{"object"}, which
  is an object that inherits from the \code{"parameter"} class. The
  \code{make_parameter} methods are defined for each inherited class and
  basically do two things. First, they do something like this\cr
  \code{object@x[object@free] <- par[object@select]}\cr
  to fill the free elements of \code{x} with corresponding values from
  \code{"par"}. Then, these methods often check whether the parameters
  are collectively admissable under the assumptions of the factor analysis
  model. If not, the \code{make_parameter} method must set the 
  \code{invalid} slot to some number that is greater than \eqn{-1.0} and not
  \eqn{0.0}. It is preferable if larger values somehow indicate more flagrant
  inadmissability, since this number is used as a return value in the 
  lexical optimization process (see \code{\link[rgenoud]{genoud}}). Finally,
  the \code{make_parameter} should return \code{object}.
}
\references{
  Howe, W.G. (1955) \emph{Some Contributions to Factor Analysis}. Dissertation published
    as ORNL-1919 by Oak Ridge National Laboratory in Tennessee.

  \enc{Reiersøl}{Reiersol}, O. (1950) On the Identifiability of Parameters in Thurstone's 
     Multiple Factor Analysis. \emph{Psychometrika}, \bold{15}, 121--149.
}

\author{ Ben Goodrich }
\examples{
showClass("parameter")
showClass("parameter.cormat")
showClass("parameter.scale")
showClass("parameter.coef")
showClass("parameter.coef.nl")
showClass("parameter.coef.SEFA")
showClass("parameter.coef.SEFA.nl")
showMethods("make_parameter")
}
\keyword{classes}
\keyword{methods}
