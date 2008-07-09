\name{make_manifest}
\docType{genericFunction}
\docType{methods}
\alias{make_manifest}
\alias{make_manifest-methods}
\alias{make_manifest,missing,missing,list-method}
\alias{make_manifest,missing,missing,hetcor-method}
\alias{make_manifest,missing,missing,matrix-method}
\alias{make_manifest,missing,missing,CovMcd-method}
\alias{make_manifest,data.frame,missing,missing-method}
\alias{make_manifest,missing,data.frame,missing-method}
\alias{make_manifest,missing,matrix,missing-method}
\alias{make_manifest,matrix,missing,missing-method}
\alias{make_manifest,formula,data.frame,missing-method}
\title{ Make an object that inherits from class "manifest" }
\description{
  This function is intended for users and sets up the left-hand side of the 
  factor analysis model and is a prerequisite for calling 
  \code{\link{make_restrictions}} and \code{\link{Factanal}}.

  Although it is possible to simply estimate and use the unbiased sample 
  covariance matrix, there are many other ways to estimate a covariance 
  that can be superior, particularly when the traditional maximum likelihood 
  discrepancy function is not chosen in the call to \code{\link{make_restrictions}}.

  In technical terms, \code{make_manifest} is the constructor for objects of
  \code{\link{manifest-class}}, which houses the sample covariance estimate
  and some ancillary information in its slots. The three arguments in the signature
  of the S4 generic function are:\cr
  \code{x}, \code{data}, and \code{covmat}
}
\usage{
\S4method{make_manifest}{missing,missing,list}(covmat, n.obs = NA_integer_, shrink = FALSE)
\S4method{make_manifest}{missing,missing,hetcor}(covmat, shrink = FALSE)
\S4method{make_manifest}{missing,missing,matrix}(covmat, n.obs = NA_integer_, shrink = FALSE, sds = NULL)
\S4method{make_manifest}{missing,missing,CovMcd}(covmat)

# Use the methods above when only the covariance matrix is available
# Use the methods below when the raw data are available (preferable)

\S4method{make_manifest}{data.frame,missing,missing}(x, subset, shrink = FALSE, 
bootstrap = 0, how = "default", seed = 12345, wt = NULL, ...)

\S4method{make_manifest}{missing,data.frame,missing}(data, subset, shrink = FALSE,
bootstrap = 0, how = "default", seed = 12345, wt = NULL, ...)

\S4method{make_manifest}{missing,matrix,missing}(data, subset, shrink = FALSE,
bootstrap = 0, how = "default", seed = 12345, wt = NULL, ...)

\S4method{make_manifest}{matrix,missing,missing}(x, subset, shrink = FALSE,
bootstrap = 0, how = "default", seed = 12345, wt = NULL, ...)

\S4method{make_manifest}{formula,data.frame,missing}(x, data, subset, shrink = FALSE, na.action = "na.pass", 
bootstrap = 0, how = "default", seed = 12345, wt = NULL, ...)
}
\arguments{
\item{x}{a formula, data.frame, nonsquare matrix of observations by variables, or missing.
  If a formula, then \code{data} must be a data.frame and the formula should not have a
  response. If a data.frame or a matrix of data, then all its columns are used.}
\item{data}{ a data.frame, nonsquare matrix of observations by variables, or missing. If a
  data.frame and \code{formula} is not specified, then all its columns are used and similarly
  if it is a matrix of data.}
\item{covmat}{A covariance matrix, a list, an object of \code{\link[rrcov]{CovMcd-class}}, an
  object of S4 class \code{"hetcor"} from the \pkg{polycor} package, or missing. If a list, it
  must contain an element named \code{"cov"} and may contain the following named elements:
  \describe{
     \item{n.obs}{the number of observations used in calculating the \code{"cov"} element}
     \item{W}{a positive definite matrix to be used as a weight matrix in the ADF discrepancy
       function. However, the \code{\link{make_restrictions-methods}} can calculate various
       weight matrices if the raw data are passed to \code{make_manifest}, so this mechanism
       should only be used if those options are inadequate}
     \item{sds}{a numeric vector of standard deviations to be used if \code{"cov"} is really
       a correlation matrix}
  }}
\item{n.obs}{ The number of observations, which is used if \code{covmat} is a covariance
  matrix or if \code{covmat} is a list with no element named \code{n.obs}. It is possible 
  to obtain maximum likelihood estimates without knowing the number of observations but nothing
  else}
\item{shrink}{A logical indicating whether to use a \dQuote{shrinkage} estimator of
     the covariance matrix. If \code{TRUE}, then the \dQuote{minimax shrinkage} estimator 
     discussed in theorem 3.1 of Dey and Srinivasan (1985) is invoked on the sample covariance 
     matrix as calculated according to the other arguments. In some circumstances, \code{shrink} 
     is inappropriate and ignored with a warning}
\item{sds}{Either \code{NULL} or a numeric vector that contains the standard deviations of 
     the manifest variables, which is used when \code{covmat} is a correlation matrix}
\item{subset}{ A specification of the cases to be used}
\item{bootstrap}{A nonnegative integer (defaulting to zero) indicating how many bootstraps 
     to do when estimating the uncertainty of the sample covariance estimates.}
\item{how}{A character string indicating how the covariance matrix should be estimated;
     see the Details section}
\item{seed}{A vector of length at most one to be used as the random number generator
     seed if \code{how = "mcd"} or \code{bootstrap > 0}. If \code{NULL}, then the
     current seed is used. This argument defaults to \code{12345}.}
\item{wt}{An \strong{optional} numeric vector of weights that is the same length as 
     the number of observations that indicates the weight for each observation when
     \code{x} is specified. By default, the observations are weighted equally.
     The \code{wt} argument can be used in two ways. First, it is passed to the
     the corresponding argument of \code{\link{cov.wt}} if appropriate (see below).
     Second, it is passed to the \code{prob} argument of \code{\link{sample}} when
     \code{bootstrap > 0}.}
\item{na.action}{ The \code{\link{na.action}} to be used if \code{x} is a formula.}
\item{\dots}{ Further arguments that are passed to downstream functions when \code{covmat}
  is unspecified, implying that the raw data are being used to estimate the sample covariance.}
}
\details{
The rules governing the calculation of the sample covariance matrix are as follows and
primarily depend on whether any of the manifest variables are ordered factors. First,
consider the case where all manifest variables are numeric. If any of these manifest
variables contain missing values, then the covariance matrix is estimated via maximum
likelihood under multivariate normality assumptions but requires the suggested \pkg{mvnmle}
package. Otherwise, the \code{how} argument dictates how the covariance matrix is estimated. 
There is much to be said in favor the Minimum Covariance Determinant (\code{\link[rrcov]{CovMcd}}) 
estimator (see Pison et. al. 2003) and it is used as the default when there are no missing
data, although it can subtly affect the sampling distributions of estimates that subsequently 
derived from it. The same could probably be said for the shrinkage estimators 
(either via \code{how = "lambda"} or \code{shrink = TRUE}). The Dey and Srinivasan (1985) 
shrinkage estimator preserves the eigenvectors of the preliminarily-calculated covariance 
matrix but deterministically compresses the eigenvalues. The \code{\link[corpcor]{cov.shrink}}
estimator in the \pkg{corpcor} package is based on the idea that the amount shrinkage should be
proportional to the variance of the covariance estimates. Use \code{how = "mle"} or 
\code{how = "unbiased"} to obtain either the maximum likelihood or unbiased sample covariance 
estimator, the latter of which is the one used in virtually all factor applications whether 
appropriate or not.

Next, consider the case where at least one manifest variable is an ordered factor. If
\code{how = "ranks"}, Spearman correlations are estimated from the integer codes
underlying the ordered factors. This mechanism is recommended only if there are at least
five levels of each ordered factor and no missing data. In that case, one would presumably 
want to specify \code{method = "ADF"} in the subsequent call to \code{\link{make_restrictions}}).
If \code{how != "ranks"} \strong{all} pairwise correlations are estimated under bivariate 
normality assumptions via \code{\link[polycor]{hetcor}} in the suggested \pkg{polycor} package, 
which will allow pairwise-deletion when there are missing data. If \code{how != "ranks"} and
\code{bootstrap > 0} (recommended), then there must not be any missing data because
the bootstrapping utilizes fast Spearman correlations and then tries to correct the bias
by rescaling the bootstrapped means to equal to point estimates calculated with the call
to \code{\link[polycor]{hetcor}}.

In general, bootstrapping is good for estimating the uncertainty of the estimated sample
covariances and this uncertainty estimate is needed for the ADF discrepancy function and
its special cases. In some cases, bootstrapping is the only way to obtain such an uncertainty
estimate.
}
\value{
  An object that inherits from \code{\link{manifest-class}}.
}
\author{Ben Goodrich}
\seealso{ \code{\link{Factanal}}, \code{\link{make_restrictions}}, \code{\link{manifest-class}},  
\code{\link[robustbase]{covMcd}}, \code{\link{cov.wt}}, \code{\link[polycor]{hetcor}}, 
\code{\link[mvnmle]{mlest}}, \code{\link[corpcor]{cov.shrink}}, and \code{\link{cov}}. }
\references{
  Dey, D. K. and Srinivasan K. (1985)
  Estimation of a covariance matrix under Stein's loss.
  \emph{The Annals of Statistics}, \bold{13}, 1581--1591.

  Pison, G., Rousseeuw, P.J., Filzmoser, P. and Croux, C. (2003)
  Robust factor analysis. \emph{Journal of Multivariate Analysis},
  \bold{84}, 145--172.
}
\examples{
man <- make_manifest(covmat = Harman23.cor)
show(man)      # some basic info
if(require(nFactors)) screeplot(man) # advanced Scree plot
cormat(man)    # sample correlation matrix
}
\keyword{ multivariate }
\keyword{ methods }