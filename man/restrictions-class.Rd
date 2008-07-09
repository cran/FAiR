\name{restrictions-class}
\docType{class}
\alias{restrictions-class}
\alias{restrictions.independent-class}
\alias{restrictions.factanal-class}
\alias{restrictions.orthonormal-class}
\alias{restrictions.1storder-class}
\alias{restrictions.1storder.EFA-class}
\alias{restrictions.general-class}
\alias{restrictions.2ndorder-class}
\alias{show,restrictions-method}
\title{Class "restrictions"}
\description{
  It is not necessary to understand this help page if one merely
  wants to estimate a factor analysis model. This help page is 
  intended for those who want to modify or extend FAiR or otherwise
  want some idea of how \pkg{FAiR} works \dQuote{behind the scenes}.

  This class contains information in its slots about what restrictions 
  are placed on a factor analysis model, which defines the model to
  be estimated.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("restrictions", ...)}.
However, it is not recommended to do so in normal usage because 
\code{\link{make_restrictions}} provides the constructor for users.
}
\section{Slots}{
    The \code{"restrictions"} class is virtual has the following slots:
 \describe{
    \item{\code{factors}:}{An integer vector of length two indicating the number 
      of factors at level one and level two of the model, must be nonnegative}
    \item{\code{nvars}:}{An integer indicating the number of free parameters in the 
      model, which corresponds to the \code{nvars} argument to 
      \code{\link[rgenoud]{genoud}}.}
    \item{\code{dof}:}{An integer indicating the number of degrees of freedom
      for the factor analysis model}
    \item{\code{Domains}:}{A numeric matrix with \code{nvars} rows and exactly two 
      columns indicating the lower and upper bounds for each free parameter, 
      which corresponds to the \code{Domains} argument to \code{\link[rgenoud]{genoud}}. }
    \item{\code{model}:}{A character string indicating whether a \code{"SEFA"}, 
      \code{"EFA"}, or \code{"CFA"} model is to be estimated.}
    \item{\code{discrepancy}:}{A character string indicating which discrepancy function 
       is used to estimate the model, see the options for the same argument to
       \code{\link{make_manifest}}}
    \item{\code{free}:}{A logical vector that indicates which of the parameters
      in the various matrices are free and which are considered fixed, excluding
      parameters that are \dQuote{obviously} fixed.}
  }
  An object of class \code{"restrictions.independent"} inherits from 
  \code{"restrictions"} and is used to define a \dQuote{null} model where there are
  no factors (hence \code{factors = c(0L,0L)} ) and the only parameters to estimate
  are the manifest standard deviatiations. It has two additional slots:
  \describe{
    \item{\code{scale}:}{an object of \code{\link{parameter.scale-class}}}
    \item{\code{criteria}:}{A list with one function, the discrepancy function that will
      be minimized to estimate the model and is that named by the \code{discrepancy} slot}
   }
   An object of class \code{"restrictions.factanal"} inherits from \code{"restrictions"}
   but has the same slots. This class is used for exploratory factor analysis via the
   same restrictions that are used in \code{\link{factanal}}. An object of class 
   \code{"restrictions.orthonormal"} inherits from \code{"restrictions"} and also implies 
   exploratory factor analysis but instead imposes the restrictions that the upper 
   triangle of the coefficient matrix has all zeros. It has the following additional slots:
  \describe{
    \item{\code{beta}:}{an object of \code{\link{parameter.coef-class}} where the \code{x} 
       slot is a primary pattern matrix with \eqn{0.0} in its upper triangle and \code{NA} 
       elsewhere}
    \item{\code{scale}:}{an object of \code{\link{parameter.scale-class}}}
    \item{\code{criteria}:}{A list with one function, the discrepancy function that will
      be minimized to estimate the model and is that named by the \code{discrepancy} slot}
  }
  An object of class \code{"restrictions.1storder"} inherits from \code{"restrictions"}
  and is used for factor analysis models with multiple correlated factors. It has the
  following additional slots:
  \describe{
    \item{\code{Phi}:}{an object of \code{\link{parameter.cormat-class}}}
    \item{\code{beta}:}{an object that inherits from \code{parameter.coef-class}. If
      \code{model = "CFA"} it should be of \code{\link{parameter.coef-class}} or
      \code{\link{parameter.coef.nl-class}}, depending on whether there are any nonlinear
      exact restrictions to impose. If \code{model = "SEFA"}, it should be of
      \code{\link{parameter.coef.SEFA-class}} or \code{\link{parameter.coef.SEFA.nl-class}},
      again depending on whether there are any nonlinear exact restrictions to impose.
      Unlike the case for \code{"restrictions.orthonormal"}, its \code{x} slot can have 
      any combination of free and fixed cells.}
    \item{\code{coef}:}{an object of \code{\link{parameter.scale-class}}}
    \item{criteria}{A list of functions to be evaluated as lexical criteria during the
      lexical optimization process. Hence, this list may contain more than one function
      but the function named in the \code{discrepancy} slot must be the last element.
      Preceeding list elements should be functions that return \eqn{-1.0} if the constraint
      is satisfied and some number greater than \eqn{-1.0} if the constraint is not satisfied}
  }
  An object of class \code{"restrictions.1storder.EFA"} inherits from 
  \code{"restrictions.1storder"} and is used for exploratory factor analysis after the
  preliminary factors have been transformed. It has the following additional slots:
  \describe{
    \item{\code{Lambda}:}{A numeric matrix with zeros in its upper triangle containing the 
      preliminary factor loadings, in other words the \code{x} slot of the \code{beta}
      slot of an object of class \code{"restrictions.orthonormal"}}
    \item{orthogonal}{A logical indicating whether the transformation was orthogonal or
      oblique. \code{\link{Rotate}} only supports oblique transformations, but it is 
      possible to \dQuote{import} an orthogonal transformation from the \pkg{GPArotation}
      package using the \code{\link{GPA2FA}} function}
    \item{Tmat}{A square numeric matrix whose order is equal to the number of factors that
      was found during the factor transformation process. Note that this matrix should
      have unit-length \emph{columns}, whereas many textbooks instead specify unit-length
      \emph{rows}}
    \item{Tcriteria}{A list, usually containing functions that were used as lexical criteria
      during the lexical optimization process (see \code{\link{Rotate}}) If the 
      \pkg{GPArotation} package was used to find the transformation, then this list must
      be a single character string naming the analytic criterion used}
  }
  An object of class \code{"restrictions.general"} inherits from 
  \code{"restrictions.1storder"} and is used for SEFA and CFA with a single second-order 
  factor that explains the systematic correlation among the first-order factors. It has one 
  additional slot, 
  \describe{
     \item{\code{Delta}:}{an object of \code{\link{parameter.coef-class}} or 
       \code{\link{parameter.coef.nl-class}}, depending on whether there are any
       nonlinear exact restrictions to impose. Its \code{x} slot has only one column 
       but it can have any combination of free and fixed cells}
  }
  An object of class \code{"restrictions.2ndorder"} inherits from \code{"restrictions.general"}
  and is used for SEFA and CFA models with multiple second-order factors to explain the 
  systematic correlation among the first-order factors. Its \code{Delta} slot now must
  only inherit from \code{\link{parameter.coef-class}} and hence can also be 
  \code{\link{parameter.coef.nl-class}} if \code{model = "CFA"} or, if 
  \code{model = "SEFA"} can be of \code{\link{parameter.coef.SEFA-class}} or 
  \code{\link{parameter.coef.SEFA.nl-class}}, again depending on whether there are any
  nonlinear exact restrictions to impose at level two. The \code{"restrictions.2ndorder"}
  class has one additional slot,
  \describe{
    \item{\code{Xi}:}{an object of \code{parameter.cormat-class}}
  }
}
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "restrictions")}, prints the object
       in a relatively easy-to-digest format and there are tailored methods for almost
       all of the classes that inherit from the \code{"restrictions"} class}
  }
  The following methods are listed here but documented elsewhere
  \describe{
    \item{restrictions2model}{the \code{\link{restrictions2model}} methods take a grand
      parameter vector and \dQuote{fill} in the free cells of the appropriate matrices
      with a \code{"restrictions"} object}
    \item{fitS4}{the \code{\link{fitS4}} methods evaluate all the lexical criteria
      given a \dQuote{filled} object that inherits from class \code{"restrictions"}}
    \item{bfgs\_fitS4}{the \code{\link{bfgs_fitS4}} methods return the discrepancy function
      given a \dQuote{filled} object that inherits from class \code{"restrictions"}}
    \item{gr\_fitS4}{the \code{\link{gr_fitS4}} methods return the gradient of the discrepancy
      function given a \dQuote{filled} object that inherits from class \code{"restrictions"}}
    \item{bfgs\_fitS4}{the \code{\link{bfgs_fitS4}} methods return a list that is then passed
      to the \code{\link{bfgs_fitS4}} and \code{\link{gr_fitS4}} methods; see the 
      \code{BFGShelp} argument to \code{\link[rgenoud]{genoud}}}
    \item{create\_start}{the \code{\link{create_start}} methods return a matrix of starting
      values for the genetic algorithm; see also \code{\link[rgenoud]{genoud}}}
    \item{create\_FAobject}{the \code{\link{create_FAobject}} methods return an object that
      inherits from \code{\link{FA-class}} once the optimal estimates have been found}
    \item{restrictions2draws}{the \code{\link{restrictions2draws}} methods return a list
      of arrays of simulated draws from the sampling distribution of the estimated free
      parameters}
    \item{restrictions2RAM}{the \code{\link{restrictions2RAM}} methods convert an object
      that inherits from class \code{"restrictions"} from that parameterization used in
      the \pkg{FAiR} package to the reticular action model (RAM) parameterization used in 
      the \pkg{sem} package}
    \item{loadings}{the \code{\link{loadings}} methods extact the factor loadings}
    \item{coef}{the \code{\link{coef,restrictions-method}}s also extract the factor loadings}
    \item{cormat}{the \code{\link{cormat}} methods extract the correlation matrices among
      factors}
    \item{uniquenesses}{the \code{\link{uniquenesses}} methods extract the uniquenesses}
    \item{fitted}{the \code{\link{fitted,restrictions-method}}s produce a correlation or
      covariance matrix as a function of the model parameters}
    \item{df.residual}{the \code{\link{df.residual,restrictions-method}}s extract the 
      degrees of freedom for the factor analysis model}
  }
}
\details{
  These classes are fundamental in \pkg{FAiR} and fulfill two important roles. First, they
  serve as a vessel that contains the parameter matrices and other necessary information to estimate
  different factor analysis models via \code{\link{Factanal}}. During the genetic optimization
  proccess, the \code{\link{restrictions2model}} method repeatedly calls 
  \code{\link{make_parameter-methods}} to \dQuote{fill} the \dQuote{empty} \code{x} slots
  of the matrices within the \code{"restrictions"} object that inherit from 
  \code{\link{parameter-class}}. Eventually, the optimal parameters satisfying the 
  specified constraints are found by the genetic algorithm and the \code{\link{create_FAobject}} 
  creates an object that inherits from class \code{\link{FA-class}} on the basis of this 
  \dQuote{filled} \code{"restrictions"} object. 

  The second purpose of a \code{"restrictions"} object is to remind the user of exactly what 
  restrictions were imposed during the estimation process, which is why a tailored 
  \code{\link{show}} method is essential. 

  Here is a table of the class hierarchy:
\tabular{lll}{
  \strong{Class} \tab \strong{Parent Class} \tab \strong{Corresponding Model}\cr
  \code{"restrictions"}  \tab none \tab none (virtual class)\cr
  \code{"restrictions.independent"} \tab \code{"retrictions"} \tab zero factor null model\cr
  \code{"restrictions.factanal"} \tab \code{"restrictions"} \tab EFA via factanal algorithm\cr
  \code{"restrictions.orthonormal"} \tab \code{"restrictions"} \tab EFA via a different algorithm\cr
  \code{"restrictions.1storder"} \tab \code{"restrictions"} \tab Models with correlated 
    factors\cr
  \code{"restrictions.1storder.EFA"} \tab \code{"restrictions.1storder"} \tab 
    EFA after factor transformation\cr
  \code{"restrictions.general"} \tab \code{"restrictions.1storder"} \tab 
    SEFA and CFA\cr
  \code{"restrictions.2ndorder"} \tab \code{"restrictions.general"} \tab 
    SEFA and CFA\cr
}
}
\author{Ben Goodrich}
\seealso{\code{\link{make_restrictions}}}
\examples{
showClass("restrictions")
showClass("restrictions.independent")
showClass("restrictions.orthonormal")
showClass("restrictions.1storder")
showClass("restrictions.1storder.EFA")
showClass("restrictions.general")
showClass("restrictions.2ndorder")
}
\keyword{ classes }
