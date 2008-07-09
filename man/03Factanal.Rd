\name{Factanal}
\alias{Factanal}
\title{ Estimate Factor Analysis Models }
\description{
This function is intended for users and estimates a factor analysis model that has been 
set up previously with a call to \code{\link{make_manifest}} and a call to 
\code{\link{make_restrictions}}.
}
\usage{
Factanal(manifest, restrictions, scores = "none", seeds = 12345, 
lower = sqrt(.Machine$double.eps), analytic = TRUE, reject = TRUE, 
NelderMead = TRUE, impatient = FALSE, ...)
}
\arguments{
  \item{manifest}{ An object that inherits from \code{\link{manifest-class}} and is
    typically produced by \code{\link{make_manifest}}.}
  \item{restrictions}{ An object that inherits from \code{\link{restrictions-class}} and
    is typically produced by \code{\link{make_restrictions}}.}
  \item{scores}{ Type of factor scores to produce, if any. The default is \code{"none"}. 
     Other valid choices (which can be partially matched) are \code{"regression"},
     \code{"Bartlett"}, \code{"Thurstone"}, \code{"Ledermann"},
     \code{"Anderson-Rubin"}, \code{"McDonald"},\cr\code{"Krinjen"},
     \code{"Takeuchi"}, and \code{"Harman"}. See Beauducel (2007) for 
     formulae for these factor scores as well as proofs that all but
     \code{"regression"} and \code{"Harman"} produce the same
     correlation matrix. }
  \item{seeds}{ A vector of length one or two to be used as the random
     number generator seeds corresponding to the \code{unif.seed} and
     \code{int.seed} arguments to \code{\link[rgenoud]{genoud}} respectively.
     If \code{seeds} is a single number, this seed is used for both 
     \code{unif.seed} and \code{int.seed}. These seeds override the defaults
     for \code{\link[rgenoud]{genoud}} and make it easier to replicate
     an analysis exactly. If \code{NULL}, the default arguments for 
     \code{unif.seed} and \code{int.seed} as specified in \code{\link[rgenoud]{genoud}} 
     are used. \code{NULL} should be used in simulations or else they will be
     horribly wrong. }
  \item{lower}{ A lower bound. In exploratory factor analysis, \code{lower} is the
     minimum uniqueness and corresponds to the \code{'lower'} element of the list 
     specified for \code{control} in \code{\link{factanal}}. Otherwise, \code{lower}
     is the lower bound used for singular values when checking for positive-definiteness
     and ranks of matrices. If the unlikely event that you get errors referencing
     positive definiteness, try increasing the value of \code{lower} slightly. }
  \item{analytic}{ A logical (default to \code{TRUE}) indicating whether analytic gradients
     should be used as much as possible. If \code{FALSE}, then numeric gradients will be 
     calculated, which are slower and slightly less accurate but are necessary in some 
     situations and useful for debugging analytic gradients.}
  \item{reject}{Logical indicating whether to reject starting values that fail the 
    constraints required by the model; see \code{\link{create_start}}}
  \item{NelderMead}{Logical indicating whether to call \code{\link{optim}} with
    \code{method = "Nelder-Mead"} when the genetic algorithm has finished to further
    polish the solution. This option is not relevant or necessary for exploratory factor
    analysis models.} 
  \item{impatient}{Logical that defaults to \code{FALSE}. If \code{restrictions} is of
    \code{\link{restrictions.factanal-class}}, setting it to \code{TRUE} will cause
    \code{\link{factanal}} to be used for optimization instead of \code{\link[rgenoud]{genoud}}. 
    In all other situations, setting it to \code{TRUE} will use \code{\link{factanal}} to
    to generate initial communality estimates instead of the slower default mechanism.}
  \item{\dots}{ Further arguments that are passed to \code{\link[rgenoud]{genoud}}.
     The following arguments to \code{\link[rgenoud]{genoud}} are hard-coded and
     cannot be changed because they are logically required by the factor analyis
     estimator:
     \tabular{lll}{
       \strong{argument} \tab \strong{value} \tab \strong{why?}\cr
       \code{nvars} \tab \code{restrictions@nvars} \tab \cr
       \code{max} \tab \code{FALSE} \tab minimizing the objective\cr
       \code{hessian} \tab \code{FALSE} \tab we roll our own\cr
       \code{lexical} \tab \code{TRUE} (usually) \tab for restricted optimization\cr
       \code{Domains} \tab \code{restrictions@Domains} \tab \cr
       \code{data.type.int} \tab \code{FALSE} \tab parameters are doubles\cr
       \code{fn} \tab wrapper around \code{\link{fitS4}} \tab \cr
       \code{BFGSfn} \tab wrapper around \code{\link{bfgs_fitS4}} \tab \cr
       \code{BFGShelp} \tab wrapper around \code{\link{bfgs_helpS4}} \tab \cr
       \code{gr} \tab various \tab it is complicated\cr
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
       \code{BFGSburnin} \tab \eqn{-1} \tab to delay the gradient check\cr
       \code{max.generations} \tab \eqn{1000} \tab big number is often necessary\cr
       \code{project.path} \tab contains \code{"Factanal.txt"} \tab \cr
       \code{starting.values} \tab see the Details section \tab 
     }

     Other arguments to \code{\link[rgenoud]{genoud}} will take the documented
     default values unless explicitly specified. In particular, you may want to 
     change \code{wait.generations} and \code{solution.tolerance}. Also, if informative
     bounds were placed on any of the parameters in the call to \code{make_restrictions}
     it is usually preferable to specify that \code{boundary.enforcement = 2} to use
     constrained optimization in the internal calls to \code{\link[stats]{optim}}. However,
     the \code{"L-BFGS-B"} optimizer is less robust than the default \code{"BFGS"} optimizer
     and occasionally causes fatal errors, largly due to misfortune.}
}
\details{
     The call to \code{Factanal} is somewhat of a formality in the sense that most of the
     difficult decisions were already made in the call to \code{\link{make_restrictions}}
     and the call to \code{\link{make_manifest}}. The most important remaining detail is
     the specification of the values for the starting population in the genetic algorithm.

     It is not necessary to provide starting values, since there are methods for this 
     purpose; see \code{\link{create_start}}. Also, if \code{starting.values = NA}, then
     a population of starting values will be created using the typical mechanism in
     \code{\link[rgenoud]{genoud}}, namely random uniform draws from the domain of the
     parameter.

     Otherwise, if \code{reject = TRUE}, starting values that fail one or more constraints 
     are rejected and new vectors of starting values are generated until the population is 
     filled with admissable starting values. In some cases, the constraints are quite difficult 
     to satisfy by chance, and it may be more practical to specify \code{reject = FALSE} or to 
     supply starting values explicitly. If starting values are supplied, it is helpful if at 
     least one member of the genetic population satisfies all the constraints imposed on the
     model. Note the rownames of \code{restrictions@Domains}, which indicate the proper order 
     of the free parameters.

     A matrix (or vector) of starting values can be passed as \code{starting.values}.
     (Also, it is possible to pass an object of \code{\link{FA-class}} to
     \code{starting.values}, in which case the estimates from the previous call to 
     \code{Factanal} are used as the starting values.) If a matrix, it should have
     columns equal to the number of rows in \code{restrictions@Domains} in the specified
     order and one or more rows up to the number of genetic individuals in the population. 

     If \code{starting.values} is a vector, its length can be equal to the number of rows
     in \code{restrictions@Domains} in which case it is treated as a one-row matrix, or its
     length can be equal to the number of manifest variables, in which case it is passed
     to the \code{start} argument of \code{\link{create_start}} as a vector of initial
     communality estimes, thus avoiding the sometimes time-consuming process of generating
     good initial communality estimates. This process can also be accelerated by specifying
     \code{impatient = TRUE}.
}
\note{
  The underlying genetic algorithm can print a variety of output as it progresses.
  On Windows, you either have to move the scrollbar periodically to flush the output
  to the screen or disable buffering by either going to the Misc menu or by clicking 
  Control+W. The output will, by default, look something like this
  \tabular{llllll}{
         Generation \tab First \tab Second  \tab \dots \tab Last \tab Discrepancy\cr
         number \tab constraint \tab constraint \tab \tab constraint \tab function\cr

         0 \tab -1.0   \tab -1.0 \tab \dots \tab -1.0 \tab double\cr
         1 \tab -1.0   \tab -1.0 \tab \dots \tab -1.0 \tab double\cr
         \dots \tab \dots \tab \dots \tab \dots \tab \dots \tab \dots\cr
         42 \tab -1.0   \tab -1.0 \tab \dots \tab -1.0 \tab double\cr
       }
  The integer on the far left indicates the generation number. If it appears to
  skip one or more generations, that signifies that the best individual in the
  \dQuote{missing} generation was no better than the best individual in the 
  previous generation. The sequence of \eqn{-1.0} indicates that various constraints
  are being satisfied by the best individual in the generation. Some of these
  constraints are hard-coded, some are added by the choices the user makes in the call
  to \code{\link{make_restrictions}}. The curious are referred to the source code,
  but for the most part users need not worry about them provided they are \eqn{-1.0}.
  If any but the last are not \eqn{-1.0} after the first few generations, there is a
  \strong{major} problem because no individual is satisfying all the constraints.
  The last number is a double-precision number indicating the value of the discrepancy
  function. This number will decrease, sometimes painfully slowly, sometimes intermittently,
  over the generations since the discrepancy function is being minimized, subject to the
  aforementioned constraints.
}
\author{Ben Goodrich}
\value{
  An object of that inherits from \code{\link{FA-class}}.
}
\seealso{ \code{\link{make_manifest}}, \code{\link{make_restrictions}}, and 
  \code{\link{Rotate}} }
\references{
  Barthlomew, D. J.  and Knott, M. (1990) \emph{Latent Variable Analysis
    and Factor Analysis.} Second Edition, Arnold.

  Beauducel, A. (2007) In spite of indeterminancy, many common factor score
    estimates yield an identical reproduced covariance matrix.
    \emph{Psychometrika}, \bold{72}, 437--441.

  Smith, G. A. and Stanley G. (1983)
  Clocking \eqn{g}: relating intelligence and measures of timed
  performance. \emph{Intelligence}, \bold{7}, 353--368.

  Venables, W. N. and Ripley, B. D. (2002)
  \emph{Modern Applied Statistics with S.} Fourth edition.  Springer.
}
\examples{
## Example from Venables and Ripley (2002, p. 323)
## Previously from Bartholomew and Knott  (1999, p. 68--72)
## Originally from Smith and Stanley (1983)
## Replicated from example(ability.cov)

man <- make_manifest(covmat = ability.cov)

\dontrun{
## Here is the easy way to set up a SEFA model, which uses pop-up menus
res <- make_restrictions(manifest = man, factors = 2, model = "SEFA")
}

## This is the hard way to set up a restrictions object without pop-up menus
beta <- matrix(NA_real_, nrow = nrow(cormat(man)), ncol = 2)
rownames(beta) <- rownames(cormat(man))
free <- is.na(beta)
beta <- new("parameter.coef.SEFA", x = beta, free = free, num_free = sum(free))

Phi  <- diag(2)
free <- lower.tri(Phi)
Phi  <- new("parameter.cormat", x = Phi, free = free, num_free = sum(free))
res  <- make_restrictions(manifest = man, beta = beta, Phi = Phi)

# This is how to make starting values where Phi is the correlation matrix 
# among factors, beta is the matrix of coefficients, and the scales are
# the logarithm of the sample standard deviations. It is also the MLE.
starts <- c( 4.46294498156615e-01, #  Phi_{21}
             4.67036349420035e-01, # beta_{11}
             6.42220238211291e-01, # beta_{21}
             8.88564379236454e-01, # beta_{31}
             4.77779639176941e-01, # beta_{41}
            -7.13405536379741e-02, # beta_{51}
            -9.47782525342137e-08, # beta_{61}
             4.04993872375487e-01, # beta_{12}
            -1.04604290549591e-08, # beta_{22}
            -9.44950629176182e-03, # beta_{32}
             2.63078925240678e-04, # beta_{42}
             9.38038168787216e-01, # beta_{52}
             8.43618801925473e-01, # beta_{62}
             log(man@sds))         # log manifest standard deviations

sefa <- Factanal(manifest = man, restrictions = res, 
                 # NOTE: Do NOT specify any of the following tiny values in a  
                 # real research situation; it is done here solely for speed
                 starting.values = starts, pop.size = 2, max.generations = 6,
                 wait.generations = 1)
nsim <- 101 # number of simulations, also too small for real work
show(sefa)
summary(sefa, nsim = nsim)
model_comparison(sefa, nsim = nsim)

stuff <- list() # output list for various methods
stuff$model.matrix <- model.matrix(sefa) # sample correlation matrix
stuff$fitted <- fitted(sefa, reduced = TRUE) # reduced covariance matrix
stuff$residuals <- residuals(sefa) # difference between model.matrix and fitted
stuff$rstandard <- rstandard(sefa) # normalized residual matrix
stuff$weights <- weights(sefa) # (scaled) approximate weights for residuals
stuff$influence <- influence(sefa) # weights * residuals
stuff$cormat <- cormat(sefa,  matrix = "RF") # reference factor correlations
stuff$uniquenesses <- uniquenesses(sefa, standardized = FALSE) # uniquenesses
stuff$FC <- loadings(sefa, matrix = "FC") # factor contribution matrix
stuff$draws <- FA2draws(sefa, nsim = nsim) # draws from sampling distribution

if(require(nFactors)) screeplot(sefa)  # Enhanced scree plot
profile(sefa) # profile plots of non-free parameters
pairs(sefa) # Thurstone-style plot
if(require(Rgraphviz)) plot(sefa) # DAG
}
\keyword{ multivariate }
\keyword{ models }
