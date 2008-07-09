\name{read.triangular}
\alias{read.triangular}

\title{Input a Covariance or Correlation Matrix}
    
\description{
This help file and function are essentially copied from 
\code{\link[sem]{read.moments}}, which is written by John Fox.
This function makes it simpler to input covariance (preferable) or correlation
matrix to be analyzed. The matrix is input in \strong{lower-triangular} form on 
as many lines as is convenient, omitting the above-diagonal elements. The
elements on the diagonal may also optionally be omitted, in which case they
are taken to be 1. In terms of features that \pkg{FAiR} supports, it is preferable 
to input the raw data if available, followed by inputting a covariance matrix, 
and inputting a correlation matrix is a last resort.
}

\usage{
read.triangular(file = "", diag = TRUE, 
                names = paste("X", 1:n, sep = ""), ...)
}

\arguments{
  \item{file}{The (quoted) file from which to read the model specification,
    including the path to the file if it is not in the current directory. If
    \code{""} (the default), then the specification is read from the standard
    input stream, and is terminated by a blank line.}
  \item{diag}{If \code{TRUE} (the default), then the input matrix includes
    diagonal elements, which is preferable in \pkg{FAiR}.}
  \item{names}{a character vector containing the names of the variables, to
    label the rows and columns of the moment matrix.}
  \item{\dots}{further arguments passed to \code{\link{scan}}}
}

\value{
 Returns a square, symmetric covariance or correlation matrix,
 suitable for input to the \code{covmat} argument of 
\code{\link{make_manifest}}.
}

\author{John Fox, with some small changes by Ben Goodrich}

\seealso{\code{\link{make_manifest}}}

\examples{
R.DHP <- read.triangular(diag=FALSE, names=c("ROccAsp", "REdAsp", "FOccAsp", 
                "FEdAsp", "RParAsp", "RIQ", "RSES", "FSES", "FIQ", "FParAsp"),
                nlines = 9)
    .6247     
    .3269  .3669       
    .4216  .3275  .6404
    .2137  .2742  .1124  .0839
    .4105  .4043  .2903  .2598  .1839
    .3240  .4047  .3054  .2786  .0489  .2220
    .2930  .2407  .4105  .3607  .0186  .1861  .2707
    .2995  .2863  .5191  .5007  .0782  .3355  .2302  .2950
    .0760  .0702  .2784  .1988  .1147  .1021  .0931 -.0438  .2087

R.DHP
}
\keyword{ manip }

