\name{equality_restriction-class}
\docType{class}
\alias{equality_restriction-class}
\alias{show,equality_restriction-method}

\title{Class "equality\_restriction"}
\description{
  It is not necessary to understand this help page if one merely
  wants to estimate a factor analysis model. This help page is 
  intended for those who want to modify or extend FAiR or otherwise
  want some idea of how FAiR works ``behind the scenes''.

  This class defines an equality restriction for use in a semi-exploratory
  or confirmatory factor analysis model.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("equality_restriction", ...)}.
However, rarely if ever, would a user need to construct an object this way.
The \code{\link{make_restrictions}} methods do so internally.
}
\section{Slots}{
  \describe{
    \item{\code{free}:}{an integer indicating which cell of a primary pattern
      matrix is considered \dQuote{free}}
    \item{\code{fixed}:}{an integer vector indicating which cells of a primary
      pattern matrix are constrained to be equal to the \dQuote{free} cell}
    \item{\code{dims}:}{an integer vector of length two indicating the dimensions
      of the primary pattern matrix under consideration}
    \item{\code{rownames}:}{a character vector of rownames for the primary
      pattern matrix under consideration}
    \item{\code{level}:}{either \code{1L} or \code{2L} indicating the level
      of the factor analysis model under consideration}
  }
}
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "equality_restriction")}: Prints the 
      equality restriction to the screen }
  }
}
\seealso{\code{\link{restrictions.1storder-class}}, \code{\link{make_restrictions}}}
\author{Ben Goodrich}
\examples{
showClass("equality_restriction")

man  <- make_manifest(covmat = ability.cov)

## Here is how to set up an equality restriction in a SEFA model the hard way
beta <- matrix(NA_real_, nrow = nrow(cormat(man)), ncol = 2)
rownames(beta) <- rownames(cormat(man))
beta[6,2] <- Inf # "fix" beta[6,2] to Inf provisionally

## Note 11L corresponds to beta[5,2] and 12L corresponds to beta[6,2]
ER   <- new("equality_restriction", free = 11L, fixed = 12L, dims = dim(beta),
            rownames = rownames(cormat(man)), level = 1L)
free <- is.na(beta) # i.e. all but beta[6,2]
beta <- new("parameter.coef.SEFA", x = beta, free = free, num_free = sum(free),
            equalities = list(ER))

Phi  <- diag(2)
free <- lower.tri(Phi)
Phi  <- new("parameter.cormat", x = Phi, free = free, num_free = sum(free))

## Now set up a restriction to prohibit the equal coefficients from being zero
blockers <- matrix(FALSE, nrow = nrow(cormat(man)), ncol = 2)
blockers[5:6,2] <- TRUE # these two cells are hence not allowed to both be zero

res <- make_restrictions(manifest = man, beta = beta, Phi = Phi, 
                         discrepancy = "MLE", criteria = list("block_1st"),
                         methodArgs = list(blockers = blockers))

\dontrun{
## Here is the easy way to do the same thing, in the following pop-up menus
## select the options to impose equality restrictions and inequality 
## restrictions (block those coefficients from being zero in a SEFA model)
res <- make_restrictions(manifest = man, factors = 2, model = "SEFA")
}
show(res)
}
\keyword{classes}
