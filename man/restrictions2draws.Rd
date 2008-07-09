\name{restrictions2draws}
\docType{genericFunction}
\docType{methods}
\alias{restrictions2draws}
\alias{restrictions2draws-methods}
\alias{FA2draws}
\title{Draw from sampling distribution of estimates}
\description{
  It is not necessary to understand this help page if one merely
  wants to estimate a factor analysis model. This help page is 
  intended for those who want to modify or extend FAiR or otherwise
  want some idea of how \pkg{FAiR} works \dQuote{behind the scenes}.\cr

  These functions facilitate drawing from a multivariate normal distribution
  defined by the estimated free parameters and their variance-covariance
  matrix. The draw is accepted if it satisfies all constraints and
  rejected otherwise. Appropriate transformations are then applied to acceptable
  draws. This mechanism makes it possible to calculate the uncertainty
  of estimates in a general fashion. It is unlikely that a user would
  ever need to call \code{restrictions2draws} directly; it is more common to use
  \code{\link{FA2draws}}, which calls \code{restrictions2draws} internally.
}
\usage{
FA2draws(object, nsim = 1001, seed = NULL, 
         covariances = FALSE, standardized = TRUE, ...)
}
\arguments{
  \item{object}{object of \code{\link{FA-class}}}
  \item{nsim}{number of non-rejected simulations}
  \item{seed}{ either NULL or an integer that will be used in a call to set.seed before 
     simulating the free parameters. The default, \code{NULL} will not change the random 
     generator state }
  \item{covariances}{ logical indicating whether return simulations of the reproduced
      covariances, which are a function of the free parameters, rather than simulations
      of the free parameters themselves. The default, \code{TRUE}, returns the simulations
      of the free parameters, which can then be manipulated into simulations of the
      reproduced covariances }
  \item{standardized}{ logical indicating whether the simulations should be rescaled so
      that they are calibrated to a correlation matrix }
  \item{\dots}{ other arguments passed to downstream functions; not currently used}
}
\section{Methods}{
  Note \code{FA2draws} is not a S4 generic function, but it primarily exists to call
  the S4 generic function, \code{restrictions2draws}. Methods are currently only defined 
  for objects of class \code{"manifest.basic"}, which are inherited by objects of class 
  \code{"manifest.data"} and \code{"manifest.data.mcd"}. There are methods
  for each flavor of \code{\link{restrictions-class}}, except for 
  \code{"restrictions.factanal"}.\cr
  
  There are also at least four arguments that are \strong{not} part of the
  signature of \code{restrictions2draws}. The first is \code{vcov}, which is the 
  variance-covariance matrix for the free parameters and is simply extracted from
  the slot of the object of \code{\link{FA-class}}. The second through fourth arguments 
  are \code{nsim}, \code{covariances}, and \code{standardized} as documented in the
  arguments section.

  For the \code{restrictions2draws} method defined for class \code{"restrictions.orthonormal"}, 
  there are two additional arguments. The first is \code{Tmat}, which can be
  \code{NULL} but can also be an optimal transformation matrix (found by 
  \code{\link{Rotate}}). The second is \code{criteria}, which is a list of
  criteria utilized when finding the optimal transformation matrix. If both
  of these arguments are specified, the method will produce draws of the 
  transformed parameters, after generating them from a multivariate normal
  distribution. They are extracted from the object of \code{\link{FA-class}} as
  necessary.
}
\details{
  Traditional methods of estimating the uncertainty of parameter estimates are not 
  necessarily appropriate in \pkg{FAiR} because of the extensive possibilities for
  imposing inequality restrictions on the parameters during the factor extraction stage
  and / or the factor transformation stage (in the case of EFA only). Thus, traditional
  asymptotic approximations to the sampling distribution of the parameters could place
  some mass on a region of the parameter space that is inconsistent with the inequality
  restrictions. To work around this problem, these functions draw simulations of the free
  parameters from an unconstrainted multivariate normal distribution, following King (1998).
  If the draw is consistent with the inequality restrictions, it is accepted and transformed
  as necessary; otherwise it is rejected. Once \code{nsim} acceptable simulations have been
  obtained, they are returned in a list of arrays for further analysis.
}
\author{Ben Goodrich}
\seealso{ \code{\link{confint,FA-method}}, \code{\link{simulate}}, and 
          \code{\link{summary.FA-class}} }
\value{
  A named list of arrays where in each case the third dimension has extent \code{nsim}.
}
\references{
  King, G. (1998) \emph{Unifying Political Methodology: The Likelihood Theory of 
    Statistical Inference}. University of Michigan Press. 
}
\examples{
  ## See example for Factanal()
}
\keyword{methods}
