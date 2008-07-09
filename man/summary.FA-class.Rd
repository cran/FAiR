\name{summary.FA-class}
\docType{class}
\alias{summary.FA-class}
\alias{confint,summary.FA-method}
\alias{plot,summary.FA,ANY-method}
\alias{show,summary.FA-method}

\title{Class "summary.FA"}
\description{
  It is not necessary to understand this help page if one merely
  wants to estimate a factor analysis model. This help page is 
  intended for those who want to modify or extend FAiR or otherwise
  want some idea of how FAiR works ``behind the scenes''.

  This class holds the output when \code{summary} is called on an
  object of \code{\link{FA-class}}
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("summary.FA", ...)}.
However, rarely if ever, would a user need to construct an object this way.
The \code{summary} method does so internally.
}
\section{Slots}{
  \describe{
    \item{\code{restrictions}:}{Object of \code{\link{restrictions-class}}}
    \item{\code{draws}:}{List (possibly of empty) of arrays produced by a
      possible call to \code{\link{FA2draws}}}
    \item{\code{order}:}{An integer vector of the same length as the number
      of first-order factors that can be used to reorder the factors when
      printed to the screen}
    \item{\code{polarity}:}{An integer vector of the same length as the 
      number of first-order factors whose elements are either \eqn(1L) or
      \eqn{-1L} that can be used to change the direction of one or more
      factors (after they have been reordered)}
    \item{\code{conf.level}:}{A number between zero and one exclusive that
      governs the bounds of the confidence intervals}
    \item{\code{standardized}:}{A logical that indicates whether the estimates
      are standardized to the correlation metric or unstandardized on the
      covariance metric}
    \item{\code{call}:}{The call to \code{\link{Factanal}}.}
  }
}
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "summary.FA")}: Prints the 
      summary on the screen }
    \item{plot}{\code{signature(x = "summary.FA", y = "ANY")}: Produces a heatmap}
  }
}
\seealso{\code{\link{FA-class}}}
\author{Ben Goodrich}
\examples{
## See also examples for Factanal()
showClass("summary.FA")
}
\keyword{classes}
