\name{read.cefa}
\alias{read.cefa}
\alias{read.CEFA}
\title{Read Files Produced by the Comprehensive Exploratory 
       Factor Analysis (CEFA) Software Package}
\description{
  This function tries to import files produced by CEFA Version 2.0
  by Michael W. Browne, Robert Cudeck, Krishna Tateneni, and Gerhald 
  Mels so that they can then be used in R.
}
\usage{
read.cefa(file)
read.CEFA(file)
}
\arguments{
  \item{file}{Character string giving the path to the file to be imported.}
}
\details{
  \code{read.cefa} does not support importing a matrix of factor loadings,
  (Datatype 3 in CEFA) because the file does not include enough information 
  to be useful to FAiR. Instead, import the covariance matrix (Datatype 1 in
  CEFA) or the raw dataset (Datatype 2 or 4 in CEFA) and use 
  \code{\link{Factanal}} to reestimate the model. \code{read.CEFA} is just
  an alias for \code{read.cefa}. These functions have not been tested very
  much and may fail if the file being imported does not strictly adhere to
  the file format described in the CEFA documentation.
}
\value{
  If a raw dataset is imported, \code{read.cefa} returns a dataframe with
  columns as variables. Such a dataframe can be passed to the \code{data}
  argument of \code{\link{Factanal}}. If a covariance matrix is imported, 
  \code{read.cefa} returns a two-element list with the following items:
  \item{cov}{The covariance matrix}
  \item{n.obs}{Number of observations}
  This list can be passed to the \code{covmat} argument of 
  \code{\link{Factanal}}. 
}
\seealso{
  See the \pkg{foreign} library for functions to import files created by
  SPSS, SAS, etc. and \code{\link{read.table}} for a function to
  import delimited text files and \code{\link{read.fwf}} for a function
  to import fixed-width text files.
}
\references{CEFA is available for Windows (but also runs under Wine) from
            \url{http://faculty.psy.ohio-state.edu/browne/software.php}}
\author{Ben Goodrich \url{http://wiki.r-project.org/rwiki/doku.php?id=packages:cran:fair}}
\examples{
\dontrun{
  OrgComm <- read.cefa(file = file.path("Program Files", "CEFAtool" "OrgComm.inp"))
  }
}
\keyword{ manip }
