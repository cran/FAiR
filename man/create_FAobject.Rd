\name{create_FAobject}
\docType{methods}
\alias{create_FAobject}
\alias{create_FAobject-methods}
\alias{create_FAobject,restrictions,manifest.basic-method}
\alias{FA-class}
\alias{FA.EFA-class}
\alias{FA.general-class}
\alias{FA.2ndorder-class}
\alias{BIC,FA-method}
\alias{confint,FA-method}
\alias{logLik,FA-method}
\alias{plot,FA,missing-method}
\alias{profile,FA-method}
\alias{show,FA-method}
\alias{summary,FA-method}
\alias{vcov,FA-method}
\alias{predict,FA-method}
\alias{screeplot,FA-method}
\title{Class "FA" and Its Constructors}
\description{
  It is not necessary to understand this help page if one merely
  wants to estimate a factor analysis model. This help page is 
  intended for those who want to modify or extend FAiR or otherwise
  want some idea of how \pkg{FAiR} works \dQuote{behind the scenes}.

  The classes that inherit from \code{"FA"} encapsulate estimates 
  from factor analysis models. First, the constructor(s) will be
  discussed.
}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("FA", ...)}. However,
  this use of \code{new("FA", ...)} is not recommended because both 
  \code{\link{Factanal}} and \code{\link{Rotate}} provide constructors for users 
  with the help of the formal methods defined for \code{\link{create_FAobject}}.
}
\section{Slots}{
  The \code{"FA"} class is not virtual but does serve as the basis for some inherited
  classes. Its slots are:
  \describe{
    \item{loadings}{A numeric array with as many rows as their are outcome variables, 
      as many columns as there are (first-order) factors, and and five shelves. Each shelf 
      is thus a matrix and contains the estimated primary pattern (PP), primary structure (PS), 
      reference pattern (RP), reference structure (RS), and factor contribution (FC)
      matrices respectively with the dimname indicated in parentheses.}
    \item{correlations}{A numeric array with as many rows and columns as there are
      there are (first-order) factors and three shelves. Each shelf is a correlation
      matrix and contains the estimated correlations among the primary factors (PF), among 
      reference factors (RF), and between each primary factor and its corresponding
      reference factor (PR) respectively with the dimname indicated in parentheses.}
    \item{uniquenesses}{A numeric vector with as many elements as there are outcome
      variables and contains the estimated unique variances from the factor analysis model.}
    \item{scale}{A numeric vector with as many elements as there are outcome variables
      and contains the estimated standard deviations of the outcome variables in the
      factor analysis model.}
    \item{restrictions}{Object that inherits from \code{\link{restrictions-class}}}
    \item{Jacobian}{A numeric matrix that contains the derivatives of the lower triangle
      of the reproduced covariance matrix with respect to each free parameter (by column)
      at their estimates.}
    \item{vcov}{A square, numeric matrix that is the estimated variance-covariance matrix
      among the estimated parameters}
    \item{scores}{A numeric matrix, possibly with zero rows. If, in the 
      call to \code{\link{Factanal}} the user specifies that factor
      scores should be calculated, this matrix contains the factor scores
      and will have as many rows as there are observations and as many columns as 
      there are (first-order) factors.}
    \item{manifest}{Object that inherits from \code{\link{manifest-class}}}
    \item{optimization}{A list that contains what is returned by the underlying optimization 
      algorithm when called internally by \code{\link{Factanal}} and \code{\link{Rotate}}.}
    \item{call}{This slot contains the call to \code{\link{Factanal}}.}
    \item{seeds}{A numeric matrix with two columns and either one or two rows. The first
      row contains the \code{unif.seed} used by \code{\link[rgenoud]{genoud}} and the 
      \code{int.seed} in the call to \code{\link{Factanal}}. If \code{\link{Rotate}} 
      is used to transform the factors after preliminary factors have been extracted
      as part of exploratory factor analysis, this matrix has a second row containing
      \code{unif.seed} and \code{int.seed} used in the call to \code{\link{Rotate}}.}
  }
  An object of \code{"FA.EFA"} inherits from the \code{"FA"} class and has the following
  additional slots:
  \describe{
    \item{rotated}{Logical indicating whether the loadings have been rotated}
    \item{Lambda}{A numeric matrix with the preliminary factor loadings}
    \item{trans\_mats}{A numeric array with as many rows and columns as there are factors
      and three shelves. Its first shelf is called \code{"primary"} and contains the
      transformation matrix that postmultiplies \code{Lambda} to yield the rotated 
      primary pattern matrix. Its second shelf is called \code{"reference"} and contains
      the transformation matrix that postmultiplies \code{Lambda} to yield the rotated
      reference structure matrix. Its third shelf is called \code{"T"} and contains the
      matrix whose \code{\link{crossprod}} is the correlation matrix among primary 
      factors. If \code{rotated = FALSE}, then all of these transformation matrices are
      identity matrices.}
  }
  An object of \code{"FA.general"} inherits from the \code{"FA"} class and its 
  \code{restrictions} slot has an object of \code{\link{restrictions.general-class}}
  It also has the following additional slots:
  \describe{
    \item{loadings\_2nd}{A numeric array that has the same form as that in the 
      \code{loadings} slot. However, since there is only one second-order factor,
      there is no distinction among the various pattern and structure matrices. The
      factor contribution matrix in the fifth shelf is simply the square of these loadings}
    \item{uniqunesses\_2nd}{A numeric vector giving the uniquenesses at level two}
  }
  An object of \code{"FA.2ndorder"} inherits from the \code{"FA.general"} class and its
  \code{restrictions} slot has an object of \code{\link{restrictions.2ndorder-class}}.
  It also has an additional slot, 
  \describe{
     \item{\code{correlations_2nd}:}{a numeric array whose form is similar to that in the 
       \code{correlations} slot but obviously pertains to the correlations among
       second-order factors rather than first-order factors.}
  } 
}
\section{Methods}{
  The \code{create_FAobject} methods \emph{construct} an object that inherits from class
  \code{"FA"} but their signatures hinge on \code{\link{restrictions-class}} and
  (in the future) \code{\link{manifest-class}}.
  The following methods are defined for the various classes that inherit from \code{"FA"}:
  \describe{
    \item{BIC}{\code{signature(object = "FA")}}
    \item{coef}{\code{signature(object = "FA")}}
    \item{confint}{\code{signature(object = "FA")}}
    \item{logLik}{\code{signature(object = "FA")}}
    \item{plot}{\code{signature(x = "FA", y = "ANY")}}
    \item{profile}{\code{signature(fitted = "FA")}}
    \item{show}{\code{signature(object = "FA")}}
    \item{summary}{\code{signature(object = "FA")}}
    \item{vcov}{\code{signature(object = "FA")}}
  }
  In addition, the following methods are defined for classes that inherit from class
  \code{"FA"} that are documented in \code{\link{loadings}} or \code{\link{S3methodsFAiR}} 
  but are also listed here:
  \describe{
    \item{loadings}{extracts various \dQuote{loadings} matrices, use insted of \code{coef}}
    \item{cormat}{extracts various correlation matrices}
    \item{uniquenesses}{extracts uniquenesses}
    \item{pairs}{Thurstone-style \code{\link{pairs}} plot among reference structure
      correlations taken two reference factors at a time}
    \item{fitted}{the covariance or correlation matrix as reproduced by the model}
    \item{residuals}{covariance or correlation residuals}
    \item{rstandard}{covariance residuals normalized by manifest standard deviations}
    \item{weights}{(approximate) weight matrices for each correlation}
    \item{influence}{a weighted residuals matrix}
    \item{df.residual}{extract the degrees of freedom}
    \item{deviance}{the value of the discrepancy function (scaled by one less the number
      of observations)}
    \item{model.matrix}{extract that sample covariance or correlation matrix}
    \item{simulate}{simulated draws from the sampling distribution of the covariance
      or correlation matrix as reproduced by the model}
  }
  In addition, the following functions are \emph{not} S4 generics but nevertheless expect
  an object that inherits from \code{"FA"} class and will throw an error otherwise:
  \describe{
     \item{model\_comparison}{the \code{\link{model_comparison}} function produces test
       statistics and fit indices}
     \item{paired\_comparison}{the \code{\link{paired_comparison}} function tests one 
       model against another model in which it is nested}
     \item{FA2draws}{the \code{\link{FA2draws}} function is essentially a wrapper around
       the \code{\link{restrictions2draws}} generic function but is more convenient}
     \item{FA2RAM}{the \code{\link{FA2RAM}} function is essentially a wrapper around the
       \code{\link{restrictions2RAM}} generic function but is more convenient}
     \item{Rotate}{the \code{\link{Rotate}} function finds an optimal transformation of
       preliminary factors in exploratory factor analysis}
     \item{GPA2FA}{the \code{\link{GPA2FA}} function requires an object of \code{"FA.EFA"}
       class}
  }
}
\usage{
\S4method{create_FAobject}{restrictions,manifest.basic}(restrictions, manifest, opt, call, scores, lower, 
                                                        analytic)
}
\arguments{
  \item{restrictions}{object of \code{\link{restrictions-class}}}
  \item{manifest}{object of \code{\link{manifest-class}}}
  \item{opt}{The list produced by \code{\link[rgenoud]{genoud}}}
  \item{call}{the call to \code{\link{Factanal}}}
  \item{scores}{A character string indicating what kind of factor scores to
    calculate; see the same argument to \code{\link{Factanal}}}
  \item{lower}{A small numeric scalar indicating the lower bound for positive
     definiteness or minimum uniqueness; see the corresponding argument
     to \code{\link{Factanal}}.}
  \item{analytic}{A logical indicating whether analytic gradients were used.}
}
\details{
  The methods for \code{create_FAobject} are called internally right at the end of 
  \code{\link{Factanal}}. They take the result of the optimization and produce an
  object that inherits from class \code{"FA"}, which is conceptually simple, although the 
  implementation is somewhat complicated and relies on a bunch of helper functions that 
  are not exported.
}
\value{
  The methods for \code{create_FAobject} produce an object of class \code{"FA"} or that
  inherits from class \code{"FA"} as appropriate.
}
\author{Ben Goodrich}
\seealso{\code{\link{Factanal}}}
\examples{
showClass("FA")
showClass("FA.EFA")
showClass("FA.general")
showClass("FA.2ndorder")
}
\keyword{classes}

