\name{GPA2FA}
\alias{GPA2FA}
\title{Covert Result of Gradient Projection Algorithm}
\description{
  This utility function can be used in an exploratory factor analysis
  when the factors are extracted via \code{\link{Factanal}} but transformed
  using the \code{\link[GPArotation]{GPFoblq}} or 
  \code{\link[GPArotation]{GPFoblq}} functions in the suggested
  \pkg{GPArotation} package. This function simply synthesizes the results
  of each to produce an object of \code{\link{FA-class}} that will be
  recognized by the various post-estimation methods and other functions.
}
\usage{
GPA2FA(GPAobject, FAobject)
}
\arguments{
  \item{GPAobject}{the result of a call to one of the transformation functions in
    \pkg{GPArotation} (see \code{\link[GPArotation]{GPA}}) with \code{A = loadings(FAobject)}}
  \item{FAobject}{an object of \code{\link{FA.EFA-class}} produced by \code{\link{Factanal}}}
}
\details{
  In some cases, it may be preferable to use the gradient projection algorithm in the
  \pkg{GPArotation} package rather than the genetic algorithm used by 
  \code{\link{Rotate}}. The gradient projection algorithm is faster, simpler, and 
  allows for the possibility of orthogonal rotation. If the gradient projection algorithm
  is used, then this function permits a seamless transition from an object of S3 class
  \code{"GPArotation"} back to a S4 object of \code{\link{FA.EFA-class}} so that the rest
  of the post-estimation methods and functions defined in \pkg{FAiR} can be used
  on the result.

  \code{\link{Rotate}} does not do orthogonal rotation for computational and philosophical
  reasons but the oblique criteria in \pkg{GPArotation} have now been copied into \pkg{FAiR}
  so that they can be utilized by \code{\link{Rotate}}. The genetic algorithm permits
  constrained optimization, which makes it possible to avoid factor collapse and more 
  generally to optimize with respect to the intersection of a set of criteria. Also,
  \code{\link{Rotate}} allows one to use the reference structure or factor contribution
  matrix instead of the primary pattern matrix when optimizing with respect to one of
  the criteria in \pkg{GPArotation}.
}
\author{Ben Goodrich}
\seealso{\code{\link[GPArotation]{GPA}}}
\value{
  An object of \code{\link{FA.EFA-class}} with transformed loadings.
}
\examples{
  # step 1: extract factors with Factanal(), rather than factanal()
  man <- make_manifest(covmat = ability.cov)
  res <- make_restrictions(man, factors = 2, model = "EFA")
  efa <- Factanal(manifest = man, restrictions = res, impatient = TRUE)

  # steps 2 and 3: transform with quartimin() and then synthesize objects
  if( require(GPArotation)) {
    efa.GPA <- quartimin(loadings(efa))
    efa.Rotated <- GPA2FA(efa.GPA, efa)
  }
  if(!require(GPArotation)) { # steps 2 and 3 are equivalent to this
    efa.Rotated <- Rotate(efa, criteria = list("quartimin"), methodArgs = 
                          list(nfc_threshold = 0, matrix = "PP"))
  }

  # step 4: interpret with various commands, such as
  summary(efa.Rotated)
}
\keyword{ models }
\keyword{ multivariate }
