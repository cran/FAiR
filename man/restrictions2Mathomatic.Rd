\name{restrictions2Mathomatic}
\alias{restrictions2Mathomatic}
\alias{restrictions2mathomatic}
\title{Format Equations for Mathomatic}
\description{
  This function exports a system of equations that can be read by the 
  program Mathomatic by George J. Gesslein II.
}
\usage{
restrictions2Mathomatic(object, file = "")
restrictions2mathomatic(object, file = "")
}
\arguments{
  \item{object}{An object of \code{\link{restrictions-class}} or
    \code{\link{FA-class}}.}
  \item{file}{Character string giving the path to the export file. This
    argument is passed to \code{\link{cat}}, so by default the file
    is printed to the screen (and may be wrapped incorrectly).}
}
\details{
  \code{restrictions2mathomatic} is just an alias for \code{restrictions2Mathomatic},
  which \dQuote{symbolically} multiplies the correlation matrix among
  primary factors from the left and right by the primary pattern matrix and
  equates the lower triangle (exclusive of the diagonal) of this reduced
  covariance matrix to the lower triangle of the population covariance matrix.
  The resulting system of equations is exported in a format suitable for reading
  into Mathomatic, which is a lightweight program that can solve systems of
  nonlinear equations symbolically. Thus, one can use Mathomatic to verify that
  a factor analysis model is identified in the sense that (if the model holds
  exactly) the free parameters can be expressed solely in terms of the population
  covariances. See, for example, Steiger (2002) illustrations, although 
  \code{\link{Factanal}} should automatically avoid the specific problems that
  Steiger (2002) is most concerned with.\cr 

  You have to manually add equations to the system that reflect the constraints
  imposed on the model; doing so is trivial at the Mathomatic prompt. However, for a 
  variety of reasons, actually \emph{solving} the system of equations in Mathomatic
  may involve more work and frustration than one might expect. For example, if you
  have more than \eqn{n = 14} manifest variables, then you have to recompile Mathomatic
  changing \code{#define N_EQUATIONS} in \code{Mathomatic/am.h}  to some number bigger
  than \eqn{100} to hold all \eqn{0.5 * n * (n - 1)} equations implied by the factor
  analysis model, plus equations for the constraints.
}
\value{
  Nothing is returned. A file is written to \code{file}.
}
\references{Mathomatic is licensed under LGPL 2.1 and is available in source and
  binary forms from \url{http://mathomatic.orgserve.de/math/}.

  Steiger, J. H. (2002) When constraints interact: A caution about reference variables
    identification constraints, and scale dependencies in structural equation modeling.
    \emph{Psychological Methods}, \bold{7}, 210--227.
  }
\author{Ben Goodrich}
\examples{
    man <- make_manifest(covmat = ability.cov)
    beta <- matrix(NA_real_, nrow = nrow(cormat(man)), ncol = 2)
    rownames(beta) <- rownames(cormat(man))
    free <- is.na(beta)
    beta <- new("parameter.coef.SEFA", x = beta, free = free, num_free = sum(free))

    Phi  <- diag(2)
    free <- lower.tri(Phi)
    Phi  <- new("parameter.cormat", x = Phi, free = free, num_free = sum(free))
    res  <- make_restrictions(manifest = man, beta = beta, Phi = Phi,
                              discrepancy = "MLE")

    restrictions2Mathomatic(res, file = "") # file printed to screen
}
\keyword{ manip }
