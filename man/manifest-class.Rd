\name{manifest-class}
\docType{class}
\alias{manifest-class}
\alias{manifest.basic-class}
\alias{manifest.basic.userW-class}
\alias{manifest.data-class}
\alias{manifest.data.ordinal-class}
\alias{manifest.data.ranks-class}
\alias{manifest.data.mcd-class}
\alias{cormat,manifest.basic-method}
\alias{model.matrix,manifest.basic-method}
\alias{show,manifest.basic-method}
\alias{plot,manifest.basic,missing-method}
\alias{screeplot,manifest.basic-method}
\title{Class "manifest"}
\description{
  It is not necessary to understand this help page if one merely
  wants to estimate a factor analysis model. This help page is 
  intended for those who want to modify or extend FAiR or otherwise
  want some idea of how \pkg{FAiR} works \dQuote{behind the scenes}.

  This class contains information in its slots about the left-hand side
  of a factor analysis model.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("manifest", ...)}.
However, it is not recommended to do so in normal usage since 
\code{\link{make_manifest}} provides the constructor for users.
}
\section{Slots}{
    The \code{"manifest"} is virtual and has the following slots:
 \describe{
    \item{n.obs}{An integer indicating the number of observations.}
    \item{how}{A character string indicating how the covariance matrix
      among the manifest variables was calculated.}
    \item{call}{The call to \code{\link{make_manifest}}.}
  }
   An object of class \code{"manifest.basic"} inherits from the \code{"manifest"} 
   class and has five additional slots: 
  \describe{
    \item{cov}{A sample covariance  matrix among manifest variables.}
    \item{cor}{A sample correlation matrix among manifest variables.}
    \item{sds}{A numeric vector containing the square root of the
       diagonal of the sample covariance matrix.}
    \item{center}{A numeric vector containing the mean of the 
       manifest variables. }
    \item{diag}{A logical indicating whether the diagonal elements of the sample
        \dQuote{covariance} matrix are all \eqn{1.0}.}
  }
  An object of class \code{"manifest.data"} inherits from the \code{"manifest.basic"}
  class and and has four additional slots: 
  \describe{
    \item{X}{A numeric matrix of manifest variables (unstandardized).}
    \item{wt}{A numeric vector whose length is the same as the number of rows
       in \code{X} indicating the weight of each observation in calculating the 
       sample covariance matrix.}
     \item{acov}{An object of \code{\link[Matrix]{dsyMatrix-class}} from the \pkg{Matrix} 
       package representing an asymptotic estimate of the fourth-order central moments of 
       the data}
  }
  An object of class \code{"manifest.data.ordinal"} inherits from the \code{"manifest.data"}
  but has no additional slots. It is used when some of the manifest variables are ordered
  factors. Likewise for the \code{"manifest.data.ranks"} class which has no additional
  slots and is used when Spearman correlations are calculated. An object of class 
  \code{"manifest.data.mcd"} inherits from the \code{"manifest.data"} class and has one 
  additional slot, \code{CovMcd}, which is an object of \code{\link[rrcov]{CovMcd-class}}
  that is defined in the \pkg{rrcov} package.
}
\section{Methods}{
    Technically, \code{"manifest"} is in the signature of all the model-fitting generic
    functions documented in \code{\link{S4GenericsFAiR}} and the same is true for
    \code{\link{make_restrictions}} However, at the moment those methods are only defined
    for the \code{"manifest.basic"} class. The following methods are defined for the
    \code{"manifest.basic"} class, but some of them are also defined for inherited classes
    if their behavior needs to differ in small ways.
    \describe{
      \item{cormat}{\code{signature(object = "manifest.basic")}}
      \item{model.matrix}{\code{signature(object = "manifest.basic")}, which also has
         a \code{standardized = TRUE} argument that controls whether the sample
         correlation matrix or sample covariance matrix is extracted}
      \item{show}{\code{signature(object = "manifest.basic")}}
      \item{plot}{\code{signature(x = "manifest.basic", y = "ANY")}}
      \item{screeplot}{\code{signature(x = "manifest.basic")}}
    }
}
\details{
In the future, the advantages of formal S4 classes will be exploited more fully. For 
example, to include a new model in \pkg{FAiR}, the first step might be to define a 
a new class that extends \code{"manifest"} if necessary. Then, appropriate methods for 
this class would need to be written. Here is a table of the class hierarchy:
\tabular{lll}{
  \strong{Class} \tab \strong{Parent Class} \tab \strong{Comment}\cr
  \code{"manifest"}  \tab none \tab none (virtual class)\cr
  \code{"manifest.basic"} \tab \code{"manifest"} \tab Can be used for any model\cr
  \code{"manifest.basic.userW"} \tab \code{"manifest.basic"} \tab ADF with user-specified
    weight matrix\cr
  \code{"manifest.data"} \tab \code{"manifest.basic"} \tab Can be used for any model\cr
  \code{"manifest.data.ordinal"} \tab \code{"manifest.data"} \tab Used when some variables
    are ordinal\cr
  \code{"manifest.data.ranks"} \tab \code{"manifest.data"} \tab Used when variables are
    converted to their ranks\cr
  \code{"manifest.data.mcd"} \tab \code{"manifest.data"} \tab Can be used for any model
}
}
\author{Ben Goodrich}
\seealso{
\code{\link{make_manifest}} and \code{\link{S4GenericsFAiR}}
}
\examples{
showClass("manifest")
showClass("manifest.basic")
showClass("manifest.data")
showClass("manifest.data.mcd")
man <- make_manifest(covmat = Harman74.cor)
show(man)
plot(man)
str(man)
}
\keyword{ classes }
