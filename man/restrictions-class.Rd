\name{restrictions-class}
\docType{class}
\alias{restrictions-class}
\alias{restrictions.factanal-class}
\alias{restrictions.orthonormal-class}
\alias{restrictions.1storder-class}
\alias{restrictions.general-class}
\alias{restrictions.2ndorder-class}
\alias{show,restrictions-method}
\alias{show,restrictions.factanal-method}
\alias{show,restrictions.orthonormal-method}
\alias{show,restrictions.1storder-method}
\alias{show,restrictions.general-method}
\alias{show,restrictions.2ndorder-method}
\title{Class "restrictions"}
\description{
  It is not necessary to understand this help page if one merely
  wants to estimate a factor analysis model. This help page is 
  intended for those who want to modify or extend FAiR or otherwise
  want some idea of how FAiR works ``behind the scenes''.\cr

  This class contains information in its slots about what restrictions 
  are placed on a factor analysis model. As an overview, objects of class 
  \code{"restrictions.factanal"} are handled by methods that are 
  analagous to the algorithm in \code{\link[stats]{factanal}}. Objects
  of class \code{"restrictions.orthonormal"} are handled by methods that 
  estimate an exploratory factor analysis model with orthogonal factors
  and zeros in the upper triangle of the coefficient matrix. Objects of
  class \code{"restrictions.1storder"}, \code{"restrictions.general"}, and
  \code{"restrictions.2ndorder"} are handled by methods that respecitvely  
  estimate semi-exploratory or confirmatory models for a single equation model
  with correlated factors, a model with a general second-order factor
  and a model with multiple correlated second-order factors.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("restrictions", ...)}.
However, it is not recommended to do so in normal usage. \code{\link{make_restrictions}} provides a constructor that 
reacts to the responses the user gives to the pop-up menus, but it too
would rarely need to be called explicitly since it is called internally by
\code{\link{Factanal}} when the \code{restrictions} argument to \code{\link{Factanal}} is unspecified. There are several examples of manual
creation of objects of class \code{"restrictions"} in the tests.R file in the FAiR/tests subdirectory and in the examples for Factanal.Rd in the FAiR/man subdirectory.
}
\section{Slots}{
    An object of class "restrictions" has the following slots:
 \describe{
    \item{factors:}{Object of class \code{"numeric"}. A numeric vector of length two
      indicating the number of factors in the model. In the case of a one-level model
      the second element should be zero.}
    \item{nvars:}{Object of class \code{"numeric"}. The number of
      free parameter in the model, which corresponds to the \code{nvars}
      argument to \code{\link[rgenoud]{genoud}}.}
    \item{Domains:}{Object of class \code{"matrix"} A matrix with
      \code{nvars} rows and two columns indicating the lower and upper
      bounds for each free or potentially free parameter in the model,
      which corresponds to the \code{Domains} argument to 
      \code{\link[rgenoud]{genoud}}. }
    \item{model:}{Object of class \code{"character"} A character
      string indicating whether a SEFA, EFA, or CFA model is to be estimated.}
    \item{method:}{Object of class \code{"character"} A character
      string indicating whether MLE or YWLS is to be used to estimate
      the model.}
  }
   An object of class \code{"restrictions.factanal"} extends \code{"restrictions"} 
   and has one additional slot, named \code{fast}, which is a logical that indicates
   whether \code{\link[stats]{factanal}} should be called (essentially) in lieu of 
   \code{\link[rgenoud]{genoud}}.\cr
   An object of class \code{"restrictions.orthonormal"} extends \code{"restrictions"}
   and has the following additional slots:
  \describe{
    \item{Phi:}{An identity matrix with rows and columns equal to the number
      of (first order) factors that serves as the orthogonal correlation 
      matrix among the factors.}
    \item{beta:}{A list with the following named elements: \code{beta}, 
      which is a matrix with as many columns as there are outcome variables
      and as many columns as there are factors. There should be zeros in the
      upper triangle of this matrix and \code{NA_real_} elsewhere. \code{free}
      is a logical vector indicating which cells of \code{beta} are 
      \code{NA_real_}. \code{num_free} is a scalar indicating the number of
      cells in \code{beta} that are \code{NA_real_}.}
     \item{Theta2:}{A list with one named element, \code{Theta2}, which is
       a diagonal covariance matrix of uniquenesses.}
     \item{criteria:}{A list with one named element containing the function
       for which the optimal parameters are sought, such as the 
       log-likelihood function.}
  }
  An object of class \code{"restrictions.1storder"} extends \code{"restrictions"}
  and has the same additional slots as objects of class 
  \code{"restrictions.orthonormal"}. There are, however, a few differences. First, 
  \code{Phi} is a diagonal matrix with \eqn{0.5} in all of its diagonal elements. 
  Second, the \code{criteria} slot can be a list with more than one 
  function to be evaluated in the lexical optimization process. Third, the 
  \code{beta} element of the \code{beta} slot can have any pattern of 
  \code{NA_real_} and / or fixed numbers in its cells to indicate free and fixed 
  coefficients respectively. Third, the \code{beta} slot has one or two additional
  list elements. The first of these is \code{select} which is a logical
  vector indicating which elements of the grand parameter vector correspond 
  to free cells of the coefficient matrix. The second of these, which is 
  included only when semi-exploratory factor analysis is utilized, is 
  \code{fix_beta_args} and is itself a list with named elements corresponding 
  to the important but not exported function \code{FAiR:::FAiR_fix_coefficients}.\cr

  An object of class \code{"restrictions.general"} extends \code{"restrictions"} 
  and has the same additional slots as objects of class 
  \code{"restrictions.1storder"} with one addition:
  \describe{
    \item{Delta:}{A list with the following named elements: \code{Delta}, 
      which is a one-column matrix with as many rows as there are
      first-order factors factors. Free elements of this matrix should 
      be \code{NA_real_} and fixed elements should be specified as
      numbers on the \eqn{(-1,1)} interval. \code{free} is a logical vector 
      indicating which cells of \code{Delta} are  \code{NA_real_}. 
      \code{num_free} is a scalar indicating the number of
      cells in \code{Delta} that are \code{NA_real_}.}
  }
  An object of class \code{"restrictions.2ndorder"} extends \code{"restrictions"}
  and has the same additional slots as objects of class 
  \code{"restrictions.general"} with the following addition:
  \describe{
    \item{Xi:}{A diagonal matrix with as many rows and columns as there are
      second-order factors and has \eqn{0.5} in all its diagonal elements.}
   } In addition, the \code{Delta} element of the \code{Delta} slot is a
   matrix with more than one column to estimate the effect of multiple
   second-order factors with any pattern of \code{NA_real_} and fixed
   numbers to indicate free and fixed coefficients. Also, the \code{Delta} 
   slot has additional list elements, \code{select} and possibly 
   \code{fix_Delta_args}, which are analagous to the elements with the
   corresponding names in the \code{beta} slot.
}
\details{
  Essentially, an object that inherits from class \code{"restrictions"} acts as
  a vessel that contains matrices and necessary information to estimate
  the corresponding model. The \code{par} argument to the methods for
  \code{\link{fitS4-methods}}, \code{\link{bfgs_fitS4-methods}}, and
  \code{\link{gr_fitS4-methods}} fills the free elements of the 
  appropriate matrix and then the appropriate functions are evaluated in search
  of the optimal parameters that satisfy the constraints. For more details, see
  the Classes.R and S4methods.R files in the FAiR/R subdirectory.
}
\author{Ben Goodrich \url{http://wiki.r-project.org/rwiki/doku.php?id=packages:cran:fair}}
\seealso{
\code{\link{Factanal}} and \code{\link{S4GenericsFAiR}}
}
\examples{
showClass("restrictions")
showClass("restrictions.orthonormal")
showClass("restrictions.1storder")
showClass("restrictions.general")
showClass("restrictions.2ndorder")
}
\keyword{ classes }
