\name{Rotate}
\alias{Rotate}
\title{ Choose a Transformation in Exploratory Factor Analysis }
\description{
  This function should be used following \code{\link{Factanal}}
  with \code{model = "EFA"} to obtain a transformation of the 
  preliminary factors that is optimal with respect to some
  intersection of criteria.
}
\usage{
Rotate(FAobject, criteria, weights = NULL, seeds = 12345, ...)
}
\arguments{
  \item{FAobject}{ An object of class \code{"FA"} produced by
   \code{\link{Factanal}} with \code{model = "EFA"}. }
  \item{criteria}{ An \strong{optional} list whose elements are functions or
    character strings naming functions to be used in the lexical
    optimization process. See the details section. It is usually best to leave 
    this argument unspecified and allow \code{Rotate} to prompt you for the
    criteria using its pop-up menus. This argument is primarily intended 
    for use in simulations, but it is also useful for passing user-defined
    functions.\cr 

    If \code{criteria} is a list that includes character strings
    the strings should be one or more of \code{"dets"}, 
    \code{"limit_correlations"}, \code{"no_suppressors"} and 
    \code{"positive_manifold"}. The last element of the list could be a
    user-defined function or one of \code{"phi"}, \code{"varphi"},
    \code{"minimaximin"}, or \code{"LS"}.}
  \item{weights}{ This argument is only relevant if the varphi criterion
    is used for transformation and should either be NULL or a numeric vector 
    giving the weights to be used If NULL, then a dynamic weighting scheme 
    will be used for the with the varphi criterion.}
  \item{seeds}{ A vector of length one or two to be used as the random
     number generator seeds corresponding to the \code{unif.seed} and
     \code{int.seed} arguments to \code{\link[rgenoud]{genoud}} respectively.
     If \code{seeds} is a single number, this seed is used for both 
     \code{unif.seed} and \code{int.seed}. These seeds override the defaults
     for \code{\link[rgenoud]{genoud}} and make it easier to replicate
     an analysis exactly. }
  \item{\dots}{ Further arguments that are passed to \code{\link[rgenoud]{genoud}}.
     Note that several of the default arguments to \code{\link[rgenoud]{genoud}}
     are silently overridden by \code{Factanal} out of logical necissitity.
     These overridden defaults are: \code{nvars}, \code{max}, \code{hessian}, 
     \code{lexical}, \code{Domains}, \code{default.domains}, \code{data.type.int}, 
     \code{fn}, \code{BFGSfn}, \code{gr}, \code{BFGShelp}, \code{unif.seed}, 
     and \code{int.seed}. In addition, several of the default arguments to
     \code{\link[rgenoud]{genoud}} are silently overridden unless the user
     explicitly specifies them by passing the arguments through \code{dots}. 
     Such arguments include: \code{MemoryMatrix = FALSE}, \code{print.level = 1}, 
     \code{P9mix = 1}, \code{MemoryMatrix = FALSE}, \code{print.level = 1}, 
     \code{max.generations = 1000}, \code{project.path}, 
     \code{starting.values}, and \code{boundary.enforcement = 1}. The arguments to 
     \code{\link[rgenoud]{genoud}} that remain at their defaults but you may want 
     to seriously consider tweaking are \code{pop.size} and 
     \code{wait.generations}.}
}
\details{
  Factor transformation is complicated and is discussed in detail in Goodrich 
  (2008). This help page should really only be used as a reminder for what the 
  various choices are, which are normally indicated by leaving \code{criteria}
  unspecifying and responding to pop-up menus. The vignette provides a step-by-step 
  guide to navigating  the pop-up menus and formally defines the criteria to follow; 
  execute \code{vignette("FAiR")} to read it.

  The basic problem is to choose a transformation of the factors that is optimal 
  with respect to some intersection of criteria. Since the objective function is 
  vector valued, lexical optimization is used via a genetic algorithm;
  see \code{\link[rgenoud]{genoud}}.

  Currently, there are four options for the ultimate criterion in the lexical
  optimization. These include Thurstone's (1935) criterion (\code{"phi"}) and a 
  generalization of Thurstone's criterion (\code{"varphi"}) that may overcome some
  of the shortcomings with phi. In addition, there is a criterion 
  (\code{"minimaximin"}) that captures Thurstone's concept of simple structure in 
  another way, and one crition for simple loadings, called the Loading Simplicity 
  Index (\code{"LS"}) that is advocated in Lorenzo-Seva (2003). 

  There are also many options for additional criteria that precede the ultimate
  criterion in the lexical ranking. These criteria can be seen as 
  \dQuote{constraints} so that one is searching for the transformation matrix
  that yields the minimum value of the ultimate criterion subject to the 
  constraints defined by previous criteria. The only required constraint is the
  constraint that the factor intercorrelation matrix be sufficiently nonsingular.

  It is not necessary to provide starting values for the parameters, since
  there are is a mechanism to generate them automatically. But a matrix of
  starting values can be passed to through the \code{dots} to
  \code{\link[rgenoud]{genoud}}. This matrix should have rows equal to the
  \code{pop.size} argument in \code{\link[rgenoud]{genoud}} and columns equal
  the number of factors squared, which corresponds to the \code{nvars} argument
  in \code{\link[rgenoud]{genoud}}. The columns correspond to the cells of
  the transformation matrix, albeit in a somewhat convoluted order. The order
  of the parameters is row-major order with respect to the transformation
  matrix, with the exception that the parameters comprosing the last row of 
  the transformation matrix are moved to the very end. Also, in contrast to
  some texts, the transformation matrix in \code{Rotate} has unit-length
  columns, rather than unit-length rows.
}
\note{
  The underlying genetic algorithm will print a variety of output as it progresses.
  On Windows, you have to move the scrollbar periodically to flush the output
  to the screen. The output will look something like this
  \tabular{llllll}{
         0 \tab -1.0   \tab -1.0 \tab \dots \tab -1.0 \tab double\cr
         1 \tab -1.0   \tab -1.0 \tab \dots \tab -1.0 \tab double\cr
         \dots \tab \dots \tab \dots \tab \dots \tab \dots \tab \dots\cr
         437 \tab -1.0   \tab -1.0 \tab \dots \tab -1.0 \tab double\cr
       }
  The integer on the very left indicates the generation number. If it appears to
  skip one or more generations, that signifies that the best individual in the
  \dQuote{missing} generation was no better than the best individual in the 
  previous generation. The sequence of negative ones indicates that various 
  constraints specified by the user are being satisfied by the best individual in 
  the generation. The curious are referred to the source code and / or the, 
  vignette, but for the most part users need not worry about them provided they
  are \eqn{-1.0}. If any but the last are not \eqn{-1.0} after the first few
  generations, there is a problem because no individual is satisfying all the 
  constraints. The last number is a double-precision number and is typically the
  logarithm of the analytic criterion specified by the user. This number will, 
  decrease, sometimes painfully slowly, sometimes intermittently, over the 
  generations since the criterion is being minimized subject to the aforementioned
  constraints. Finally, do not be particularly concerned if there are messages
  indicating a gradient check has failed because there is no strong reason to 
  expect the gradient of the (last) criterion with respect to all the cells of the 
  transformation matrix will be particularly small 
}
\value{
  An object of formal S4 class \code{"FA"}.
}
\references{
  Goodrich, B. (2008) Analytic Transformation of Factors in FAiR. Unwritten
    manuscript.

  Lorenzo-Seva, U. (2003) A factor simplicity index.
  \emph{Psychometrika}, \bold{68}, 49--60.

  Thurstone, L. L.  (1935) \emph{The Vectors of Mind.}
    Cambridge University Press.
}
\author{Ben Goodrich \url{http://wiki.r-project.org/rwiki/doku.php?id=packages:cran:fair}}
\seealso{ \code{\link{Factanal}} }
\examples{
## See first example for Factanal()
}
\keyword{ multivariate }
\keyword{ models }
