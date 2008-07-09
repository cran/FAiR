\name{Rotate}
\alias{Rotate}
\title{ Choose a Transformation in Exploratory Factor Analysis }
\description{
  This function is intended for users and finds an optimal rotation
  of preliminary factor loadings extracted via exploratory factor
  analysis. Unlike other software, \code{Rotate} finds a rotation
  that is optimal with respect to an \strong{intersection} of
  criteria, one of which is a \dQuote{analytic} criterion and the
  rest are \dQuote{constraints}. Similar to \code{\link{make_restrictions}},
  much of the functionality of \code{Rotate} is implemented via
  pop-up menus, which is the \strong{strongly} recommended way to proceed.
  Also, the vignette details what all these options \emph{mean} in
  substantive terms; execute \code{vignette("FAiR")} to read it.
}
\usage{
Rotate(FAobject, criteria = list(), methodArgs = list(), 
       normalize = rep(1, nrow(loadings(FAobject))), seeds = 12345, 
       NelderMead = TRUE, ...)
}
\arguments{
  \item{FAobject}{ An object of \code{\link{FA.EFA-class}} produced by
   \code{\link{Factanal}}. }
  \item{criteria}{ An \strong{optional} list whose elements are functions or
    character strings naming functions to be used in the lexical
    optimization process. If unspecified, \code{Rotate} will prompt you
    with pop-up menus, which is the recommended way to proceed. See the
    Details section otherwise.}
  \item{methodArgs}{ an option list with named elements for additional
    arguments needed by some criterion functions; see the Details section}
  \item{normalize}{ a vector with the same length as the number of manifest
    variables, or a function that takes the preliminary primary pattern 
    matrix as its first argument and produces such a vector, or a character 
    string that is either \code{"kaiser"} or \code{"cureton-mulaik"}, which
    are discussed more in the Details section. The rows of the preliminary 
    primary pattern matrix are \strong{divided} by this vector before the
    optimal transformation is found. Note that row-normalization may not be 
    sensible, especially if the objective function optimizes the reference structure
    matrix or the factor contribution matrix. } 
  \item{seeds}{ A vector of length one or two to be used as the random
     number generator seeds corresponding to the \code{unif.seed} and
     \code{int.seed} arguments to \code{\link[rgenoud]{genoud}} respectively.
     If \code{seeds} is a single number, this seed is used for both 
     \code{unif.seed} and \code{int.seed}. These seeds override the defaults
     for \code{\link[rgenoud]{genoud}} and make it easier to replicate
     an analysis exactly. However, if \code{seeds = NULL}, then the default
     seeds are used, which is absolutely necessary during simulations. }
  \item{NelderMead}{Logical indicating whether to call \code{\link{optim}} with
    \code{method = "Nelder-Mead"} when the genetic algorithm has finished to further
    polish the solution.}
  \item{\dots}{ Further arguments that are passed to \code{\link[rgenoud]{genoud}}.
     Note that several of the default arguments to \code{\link[rgenoud]{genoud}}
     are silently overridden by \code{Factanal} out of logical necissitity:
     \tabular{lll}{
       \strong{argument} \tab \strong{value} \tab \strong{why?}\cr
       \code{nvars} \tab \code{FAobject@restrictions@factors[1]^2} \tab \cr
       \code{max} \tab \code{FALSE} \tab minimizing the objective\cr
       \code{hessian} \tab \code{FALSE} \tab not meaningful\cr
       \code{lexical} \tab \code{TRUE} \tab resticted optimization\cr
       \code{Domains} \tab \code{NULL} \tab \cr
       \code{default.domains} \tab \eqn{1} \tab parameters are cosines\cr
       \code{data.type.int} \tab \code{FALSE} \tab parameters are doubles\cr
       \code{fn} \tab wrapper around an internal function \tab \cr
       \code{BFGSfn} \tab wrapper around an internal function \tab \cr
       \code{gn} \tab \code{NULL} \tab analytic gradients are unknown \cr
       \code{BFGShelp} \tab wrapper around an internal function \tab \cr
       \code{unif.seed} \tab taken from \code{seeds} \tab replicability\cr
       \code{int.seed} \tab taken from \code{seeds} \tab replicability      
     }

     The following arguments to \code{\link[rgenoud]{genoud}} default to values
     that differ from those documented at \code{\link[rgenoud]{genoud}} but can
     be overridden by specifying them explicitly in the \dots :
     \tabular{lll}{
       \strong{argument} \tab \strong{value} \tab \strong{why?}\cr
       \code{boundary.enforcement} \tab \eqn{1} usually \tab \eqn{2} can cause problems\cr
       \code{MemoryMatrix} \tab \code{FALSE} \tab runs faster\cr
       \code{print.level} \tab \eqn{1} \tab output is not that helpful for \eqn{>= 2}\cr
       \code{P9mix} \tab \eqn{1} \tab to always accept the BFGS result\cr
       \code{BFGSburnin} \tab \eqn{5} \tab to have a little burnin\cr
       \code{max.generations} \tab \eqn{1000} \tab big number is often necessary\cr
       \code{project.path} \tab contains \code{"Rotate.txt"} \tab 
     }
     
     The arguments to \code{\link[rgenoud]{genoud}} that remain at their defaults but you 
     may want to consider tweaking are \code{pop.size}, \code{wait.generations}, and
     \code{solution.tolerance}.}
}
\details{
  This help page should really only be used as a reminder for what the various 
  choices are, which are normally indicated by leaving \code{criteria} and 
  \code{methodArgs} unspecified and responding to pop-up menus. The vignette provides 
  a step-by-step guide to the pop-up menus and formally defines the criteria; execute 
  \code{vignette("FAiR")} to read it.

  The basic problem is to choose a transformation of the factors that is optimal 
  with respect to some intersection of criteria. Since the objective function is 
  vector valued, lexical optimization is performed via a genetic algorithm, which is
  tantamount to constrained optimization; see \code{\link[rgenoud]{genoud}}.

  The following functions can be named as constraints but must not be the last element of
  \code{criteria}:
\tabular{lll}{
  \strong{name} \tab \strong{methodArgs} \tab \strong{reminder of what function does}\cr
  \code{"no_factor_collapse"} \tab \code{nfc_threshold} \tab to prevent factor collapse\cr
  \code{"limit_correlations"} \tab \code{lower} and \code{upper} \tab limits factor intercorrelations\cr
  \code{"positive_manifold"}  \tab \code{pm_threshold} \tab forces \dQuote{positive} manifold\cr
  \code{"ranks_rows_1st"}  \tab \code{row_ranks} \tab row-wise ordering constraints\cr
  \code{"ranks_cols_1st"}  \tab \code{col_ranks} \tab column-wise ordering constraints\cr
  \code{"indicators_1st"} \tab \code{indicators} \tab designate which is the best indicator of a factor\cr
  \code{"evRF_1st"}  \tab none \tab restrict effective variance of reference factors\cr
  \code{"evPF_1st"}  \tab none \tab restrict effective variance of primary factors\cr  
  \code{"h2_over_FC_1st"}  \tab none \tab communalities \eqn{>=} factor contributions\cr
  \code{"no_neg_suppressors_1st"}  \tab \code{FC_threshold} \tab no negative suppressors\cr
  \code{"gv_1st"}  \tab none \tab generalized variance of primary <= reference factors\cr
  \code{"distinguishability_1st"}  \tab none \tab best indicators have no negative suppressors\cr
}
 
  In fact, \code{"no_factor_collapse"} is always included and is listed above only to 
  emphasize that one must specify \code{methodArgs$nfc_threshold} to avoid seeing the
  associated pop-up menu. This restriction to avoid factor collapse makes it possible
  to utilize one of the following  \dQuote{analytic} criteria that would otherwise
  result in factor collapse much of the time. One of the following can be named as the
  \strong{last} element of \code{criteria}:
\tabular{ll}{
  \strong{name} \tab \strong{methodArgs} \cr
  \code{"phi"} \tab \code{c}\cr
  \code{"varphi"} \tab \code{weights}\cr
  \code{"LS"} \tab \code{eps}, \code{scale}, and \code{E}\cr
  \code{"minimaximin"} \tab \cr
  \code{"geomin"} \tab \code{delta} \cr
  \code{"quartimin"} \tab \cr
  \code{"target"} \tab \code{Target} \cr
  \code{"pst"} \tab \code{Target} \cr
  \code{"oblimax"} \tab \cr
  \code{"simplimax"} \tab \code{k} \cr
  \code{"bentler"} \tab \cr
  \code{"cf"} \tab \code{kappa} \cr
  \code{"infomax"} \tab \cr
  \code{"mccammon"} \tab \cr
  \code{"oblimin"} \tab \code{gam} \cr
}
  
  The first four are defined only on the reference structure matrix. The remainder are copied
  from the \pkg{GPArotation} package and have the same arguments, with the exception of
  \code{"pst"} which uses \code{NA} in the target matrix for untargeted cells rather than
  also specifying a weight matrix (which is called \code{W} in the \pkg{GPArotation} package). 
  In addition, one can specify \code{methodArgs$matrix} as one of \code{"PP"}, \code{"RS"}, or 
  \code{"FC"} to use the primary pattern, reference structure, or factor contribution matrix in 
  conjunction with the criteria from the \pkg{GPArotation} package, although these criteria
  are technically defined with respect to the primary pattern matrix.

  Row-standardization should not be necessary. Row-standardization was originally
  intended to counteract some tendencies in the transformation process that can
  now be accomplished directly through lexical (i.e. constrained) optimization.
  Nevertheless, Kaiser normalization divides each row of the preliminary primary
  pattern matrix by its length and Cureton-Mulaik normalization favors rows that
  are thought to have only one large loading after transformation. Both schemes
  are thoroughly discussed in Browne (2001), which also discusses most of the
  continuous analytic criteria available in \pkg{FAiR} with the exceptions of
  those in Thurstone (1935) and Lorenzo-Seva (2003).

  It is not necessary to provide starting values for the parameters. But a matrix of
  starting values can be passed to through the \code{dots} to \code{\link[rgenoud]{genoud}}. 
  This matrix should have rows equal to the \code{pop.size} argument in 
  \code{\link[rgenoud]{genoud}} and columns equal the number of factors squared. The columns 
  correspond to the cells of the transformation matrix in column-major order. In contrast to
  some texts, the transformation matrix in \code{Rotate} has unit-length columns, rather 
  than unit-length rows.
}
\note{
  The underlying genetic algorithm will print a variety of output as it progresses.
  On Windows, you have to move the scrollbar periodically to flush the output
  to the screen. The output will look something like this
  \tabular{llllll}{
         Generation \tab First \tab Second  \tab \dots \tab Last \tab Analytic\cr
         number \tab constraint \tab constraint \tab \tab constraint \tab criterion\cr

         0 \tab -1.0   \tab -1.0 \tab \dots \tab -1.0 \tab double\cr
         1 \tab -1.0   \tab -1.0 \tab \dots \tab -1.0 \tab double\cr
         \dots \tab \dots \tab \dots \tab \dots \tab \dots \tab \dots\cr
         42 \tab -1.0   \tab -1.0 \tab \dots \tab -1.0 \tab double\cr
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
  (logarithm of the) analytic criterion specified by the user. This number will, 
  decrease, sometimes painfully slowly, sometimes intermittently, over the 
  generations since the criterion is being minimized subject to the aforementioned
  constraints. Finally, do not be particularly concerned if there are messages
  indicating a gradient check has failed because there is no strong reason to 
  expect the gradient of the (last) criterion with respect to all the cells of the 
  transformation matrix will be particularly small 
}
\author{Ben Goodrich}
\value{
  An object of \code{\link{FA.EFA-class}}.
}
\references{
  Browne, M.W. (2001) An overview of analytic rotation in exploratory factor analysis.
    \emph{Multivariate Behavioral Research}, \bold{36}, 111--150.

  Lorenzo-Seva, U. (2003) A factor simplicity index.
  \emph{Psychometrika}, \bold{68}, 49--60.

  Smith, G. A. and Stanley G. (1983)
  Clocking \eqn{g}: relating intelligence and measures of timed
  performance. \emph{Intelligence}, \bold{7}, 353--368.

  Thurstone, L. L.  (1935) \emph{The Vectors of Mind.}
    Cambridge University Press.

  Venables, W. N. and Ripley, B. D. (2002)
  \emph{Modern Applied Statistics with S.} Fourth edition.  Springer.

}
\seealso{ \code{\link{Factanal}} }
\examples{
## Example from Venables and Ripley (2002, p. 323)
## Previously from Bartholomew and Knott  (1999, p. 68--72)
## Originally from Smith and Stanley (1983)
## Replicated from example(ability.cov)

man <- make_manifest(covmat = ability.cov)
res <- make_restrictions(man, factors = 2, model = "EFA")
efa <- Factanal(manifest = man, restrictions = res, impatient = TRUE)

show(efa); summary(efa)

# 'criteria' and 'methodArgs' would typically be left unspecified 
#  and equivalent choices would be made from the pop-up menus
efa.rotated <- Rotate(efa, criteria = list("phi"),
                      methodArgs = list(nfc_threshold = 0.25, c = 1.0))
summary(efa.rotated)
pairs(efa.rotated)
## See the example for Factanal() for more post-estimation commands
}
\keyword{ multivariate }
\keyword{ models }
