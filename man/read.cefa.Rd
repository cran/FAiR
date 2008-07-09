\name{read.cefa}
\alias{read.cefa}
\alias{read.CEFA}
\alias{write.cefa}
\alias{write.CEFA}
\title{Read and Write Files Produced by or for CEFA}
\description{
  Compatability functions between R and Comprehensive Exploratory 
  Factor Analysis (CEFA) Version 2.0 by 
  Michael W. Browne, Robert Cudeck, Krishna Tateneni, and Gerhald Mels.
}
\usage{
read.cefa(file)
read.CEFA(file)
write.cefa(file, FAobject)
write.CEFA(file, FAobject)
}
\arguments{
  \item{file}{Character string giving the path to the file to be imported, 
    possibly \code{""} to print the file to the screen.}
  \item{FAobject}{object of \code{\link{FA-class}} produced by 
     \code{\link{Factanal}} under EFA}
}
\details{
  \code{read.cefa} does not support importing a matrix of factor loadings,
  (Datatype 3 in CEFA) because the file does not include enough information 
  to be useful to \pkg{FAiR}. Instead, import the covariance matrix (Datatype 1 in
  CEFA) or better the raw dataset (Datatype 2 or 4 in CEFA) and use 
  \code{\link{Factanal}} to reestimate the model.\cr

  In contrast, the only export method that is currently supported is that for
  preliminary factor loadings in a EFA. \code{read.CEFA} is just
  an alias for \code{read.cefa} and similarly for \code{write.CEFA}.\cr

  These functions have not been tested very much.
}
\value{
  If a raw dataset is imported, \code{read.cefa} returns a dataframe with
  columns as variables. Such a dataframe can be passed to the \code{data} or
  \code{x} argument of \code{\link{make_manifest}}. If a covariance matrix is 
  imported, \code{read.cefa} returns a two-element list with the following items:
  \item{cov}{The covariance matrix}
  \item{n.obs}{Number of observations}
  This list can be passed to the \code{covmat} argument of 
  \code{\link{make_manifest}}. \code{write.cefa} does not produce anything but
  writes a file to the specified location.
}
\seealso{
  \code{\link{read.triangular}} imports a triangular covariance matrix,
  the \pkg{foreign} library has functions to import files created by
  SPSS, SAS, etc., \code{\link{read.table}} imports delimited text files and 
  \code{\link{read.fwf}} imports fixed-width text files.
}
\references{CEFA is available for Windows (but also runs under Wine in Linux) from
            \url{http://faculty.psy.ohio-state.edu/browne/software.php}}
\author{Ben Goodrich}
\examples{
  \dontrun{
    OrgComm <- read.cefa(file = file.path("Program Files", "CEFAtool" 
                                          "OrgComm.inp"))
  }
}
\keyword{ manip }
