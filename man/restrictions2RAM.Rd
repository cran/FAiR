\name{restrictions2RAM}
\docType{genericFunction}
\docType{methods}
\alias{restrictions2RAM}
\alias{restrictions2RAM-methods}
\alias{FA2RAM}
\title{Convert to reticular action model}
\description{
  These functions convert to the reticular action model (RAM) format utilized by the
  \pkg{sem} package.
}
\usage{
FA2RAM(FAobject)
}
\arguments{
  \item{FAobject}{object of \code{\link{FA-class}}}
}
\section{Methods}{
  Note \code{FA2RAM} is not a S4 generic function, but it primarily exists to call
  the S4 generic function, \code{restrictions2RAM}. There are methods for each flavor
  of \code{\link{restrictions-class}}, except for \code{"restrictions.factanal"}.

  The \code{restrictions2RAM} S4 generic function takes one argument, 
  \code{restrictions}, which is an object that inherits from
  \code{\link{restrictions-class}}.
}
\section{Warning}{
  It is not possible to map exactly to the RAM specification used in the \pkg{sem}
  package for models with two levels. \code{FA2RAM} will produce a warning in such
  cases indicating that the normalizations of the latent variables differ from the
  normalizations used in \pkg{FAiR}. The reproduced covariance matrix should be
  the same (up to some numerical error) but the parameter values at the optimum
  will differ.
}
\author{Ben Goodrich}
\seealso{ See \code{\link[sem]{specify.model}} and \code{\link[sem]{sem}} in the \pkg{sem} 
  package for details on the RAM specification.
}
\value{
  A character matrix with three columns of S3 class \code{"mod"}, which is suitable for
  passing to \code{\link[sem]{sem}}.
}
\examples{
man  <- make_manifest(covmat = ability.cov)

## Here is also an example of how to set up a CFA model the hard way
beta <- matrix(NA_real_, nrow = nrow(cormat(man)), ncol = 2)
rownames(beta) <- rownames(cormat(man))
beta[2:3,1] <- beta[5:6,2] <- 0
free <- is.na(beta)
beta <- new("parameter.coef", x = beta, free = free, num_free = sum(free))

Phi  <- diag(2)
free <- lower.tri(Phi)
Phi  <- new("parameter.cormat", x = Phi, free = free, num_free = sum(free))

res <- make_restrictions(manifest = man, beta = beta, Phi = Phi, 
                         discrepancy = "MLE")
show(res)
RAM <- restrictions2RAM(res)
if(require(sem)) {
  RAM[,3] <- c(0.5, 0, 0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0, 0, 
               1, 0, 1,  0.5, 0.5, 0.5, 0.5, 0.5, 0.5)
  cfa <- sem(RAM, model.matrix(man, standardize = FALSE), man@n.obs)
  summary(cfa)
}
}
\keyword{manip}
\keyword{multivariate}
