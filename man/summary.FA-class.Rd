\name{summary.FA-class}
\docType{class}
\alias{summary.FA-class}
\alias{plot,summary.FA,ANY-method}
\alias{show,summary.FA-method}

\title{Class "summary.FA"}
\description{This class contains information produced by the \code{summary}
  method for objects of class \code{"FA"} or that inherit from class 
  \code{"FA"}.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("summary.FA", ...)}.
However, rarely if ever, would a user need to construct an object this way.
The \code{summary} methods do so internally a produce an object of class
\code{"summary.FA"})
}
\section{Slots}{
	 \describe{
    \item{\code{call}:}{Object of class \code{"language"} The call to
      \code{\link{Factanal}}.}
    \item{\code{point_estimates}:}{Object of class \code{"list"} A list
       containing the various point estimates. }
    \item{\code{zstats}:}{Object of class \code{"list"} A list containing
       the z-statistics for the various estimates. }
  }
}
\section{Methods}{
  \describe{
    \item{plot}{\code{signature(x = "summary.FA", y = "ANY")}: Does not
      currently work }
    \item{show}{\code{signature(object = "summary.FA")}: Prints the 
      summary on the screen }
	 }
}
\section{Warning }{The code producing the z-statistics is shaky at best.}
\seealso{\code{\link{FA-class}}}
\author{Ben Goodrich \url{http://wiki.r-project.org/rwiki/doku.php?id=packages:cran:fair}}
\examples{
showClass("summary.FA")
}
\keyword{classes}
