\name{Factanal}
\alias{Factanal}
\title{ Estimate Common Factor Analysis Models }
\description{
This function estimates models for semi-exploratory factor analysis (SEFA), exploratory factor analysis (EFA), and confirmatory factor analysis (CFA) using 
a genetic algorithm.
}
\usage{
Factanal(x, factors, data = NULL, covmat = NULL, n.obs = NA, subset, 
na.action, scores = "none", seeds = 12345, lower = sqrt(.Machine$double.eps), 
model = c("SEFA", "EFA", "CFA"), method = c("MLE", "YWLS"), 
restrictions, fixed, criteria = NULL, robust.covmat = FALSE, ...)
}
\arguments{
  The first several arguments (through scores) are largely similar
  to those in \code{\link{factanal}}.
  \item{x}{ A formula or a numeric matrix or an object that can be
    coerced to a numeric matrix. This argument is required if
    \code{covmat = NULL} or if \code{robust.covmat = TRUE} and is
    always recommended if the raw data are available.}
  \item{factors}{ The number (>0) of factors to be fitted, which differs
    from the argument in \code{\link{factanal}} in that 
    \code{factors} can be a numeric vector of length two to indicate
    the number of factors to extract at level one and level two of a 
    two-level semi-exploratory or confirmatory factor analysis model.
    If a single number and greater than 2 and \code{model != "EFA"}, 
    \code{Factanal} will ask whether to estimate a second level in the 
    usual case where \code{restrictions} is unspecified. }
  \item{data}{ An optional data frame (or similar: see
    \code{\link{model.frame}}), used only if \code{x} is a formula.  By
    default the variables are taken from \code{environment(formula)}. }
  \item{covmat}{ A covariance matrix, or a covariance list as returned by
    \code{\link{cov.wt}} or similar. If the covariance matrix is really a
    correlation matrix, it is not (yet) possible to (accurately) calculate 
    some measures of uncertainty, so the covariance matrix (or list) should 
    be passed  instead of the correlation matrix if available. It is always
    recommended to use \code{x} to pass the raw data to \code{Factanal},
    but \code{covmat} is required if \code{x} is unspecified.}
  \item{n.obs}{ The number of observations, used if covmat is a covariance 
    matrix. It is possible to obtain point estimates without knowing the
    number of observations, but it is not possible to calculate measures
    of uncertainty. }
  \item{subset}{ A specification of the cases to be used, if \code{x} is 
     a matrix or formula. }
  \item{na.action}{ The na.action to be used if x is used as a formula. }
   \item{scores}{ Type of scores to produce, if any. The default is \code{"none"}. 
     Other valid choices (which can be abbreviated) are \code{"regression"},
     \code{"Bartlett"}, \code{"Thurstone"}, \code{"Ledermann"},
     \code{"Anderson-Rubin"}, \code{"McDonald"}, \code{"Krinjen"},
     \code{"Takeuchi"}, and \code{"Harman"}. See Beauducel (2007) for 
     formulae for these factor scores as well as proofs that all but
     \code{"regression"} and \code{"Harman"} produce the same
     correlation matrix.}
  \item{seeds}{ A vector of length one or two to be used as the random
     number generator seeds corresponding to the \code{unif.seed} and
     \code{int.seed} arguments to \code{\link[rgenoud]{genoud}} respectively.
     If \code{seeds} is a single number, this seed is used for both 
     \code{unif.seed} and \code{int.seed}. These seeds override the defaults
     for \code{\link[rgenoud]{genoud}} and make it easier to replicate
     an analysis exactly. }
  \item{lower}{ A lower bound. In exploratory factor analysis using the fitting
     function in \code{\link[stats]{factanal}}, this argument corresponds to the 
     \code{'lower'} element of the list specified for \code{control} in 
     \code{\link{factanal}} and indicates the smallest permissable value
     for a uniqueness. Otherwise, this argument is the lower bound used for
     eigenvalues and singular values when checking for positive-definiteness
     and ranks of matrices. If you happen to get errors referencing
     positive definiteness, try increasing the value of \code{lower}. }
  \item{model}{ A character string indicating "SEFA", "EFA", or "CFA" to
     indicate whether a semi-exploratory, an exploratory, or a confirmatory 
     factor analysis model should be estimated. Defaults to "SEFA". }
  \item{method}{ A character string indicating "MLE" or "YWLS" to
     indicate how the \code{model} should be estimated. Defaults to "MLE".
     The "YWLS" option uses Yates' (1987) weighted-least squares criterion
     as opposed to most of the weighted-least squares criteria that are
     usually mentioned in the factor analysis literature. See the warning below.}
  \item{restrictions}{ An \strong{optional} object of class "restrictions". It is
     almost always best to leave this argument unspecified to allow
     \code{Factanal} to prompt for restrictions with its pop-up menus. 
     This argument is primarily intended for use in simulations. See
     \code{\link{restrictions-class}} for more information about how
     it should be specified if you do not want to be bothered by the 
     pop-up menus.}
  \item{fixed}{ An \strong{optional} matrix or list of two matrices that specifies
     the values of certain coefficients, which would be utilized most often
     in confirmatory factor analysis and is inappropriate for exploratory
     factor analysis. However, it is never necessary to explicitly supply
     this argument because \code{Factanal} will prompt you to construct
     such a matrix or matrices in confirmatory and (mixed) semi-exploratory
     factor analysis models with a correctly formatted skeleton. \cr

     If \code{fixed} is a matrix, it should have rows equal to the number of
     outcome variables and columns equal to the number of factors at level one.
     If \code{fixed} is a list of two matrices, the first element of the list 
     corresponds to the coefficient matrix at level one and the second element
     corresponds to the coefficient matrix at level two (and should have rows
     equal to the number of first-order factors and columns equal to the number
     of second-order factors). Unrestricted coefficients should be denoted with
     \code{NA} in the appropriate row and column of the appropriate matrix.
     Restricted coefficients should be denoted with the value of the restricted
     coefficient in the appropriate row and column of the appropriate matrix.}
  \item{criteria}{ An \strong{optional} list whose elements should be functions or 
     character strings that name functions to be used as criteria during
     the lexical optimization when \code{model != "EFA"}. It is almost always
     best to leave this argument unspecified to allow \code{Factanal} to prompt 
     for these criteria with pop-up menus. This argument is primarily intended for 
     use in simulations.\cr

     If \code{criteria} is a list that includes character strings
     the strings should be one or more of \code{"no_suppressors_xxx"},
     \code{"dets_xxx"}, or \code{"cohyperplanarity"} where \code{"xxx"}
     is either \code{"1st"} or \code{"2nd"} to indicate whether the
     criterion should be applied to the first or second level of the model.
     Thus, a suffix of \code{"2nd"} is only appropriate if a two-level
     model is estimated. The function implied by \code{method} is automatically
     appended to the end of this list to serve as the ultimate lexical criterion.}
  \item{robust.covmat}{A logical indicating whether a minimum covariance
     determinant estimator of the sample covariance matrix should be used. If
     \code{TRUE}, this option requires that either \pkg{robustbase} or
     \pkg{MASS} be installed; see \code{\link[robustbase]{covMcd}} and
     \code{\link[MASS]{cov.mcd}} respectively. Also, it requires that the
     raw data be passed via the \code{x} argument. Although it defaults to
     \code{FALSE}, \code{TRUE} is probably recommended in most cases.}
  \item{\dots}{ Further arguments that are passed to \code{\link[rgenoud]{genoud}}.
     Note that several of the default arguments to \code{\link[rgenoud]{genoud}}
     are silently overridden by \code{Factanal} out of logical necissitity.
     These overridden defaults are: \code{nvars}, \code{max}, \code{hessian}, 
     \code{lexical}, \code{Domains}, \code{data.type.int}, \code{fn}, \code{BFGSfn}, 
     \code{gr}, \code{BFGShelp}, \code{unif.seed}, \code{int.seed}, and --- 
     in some cases --- \code{boundary.enforcement}. In addition, several of the 
     default arguments to \code{\link[rgenoud]{genoud}} are silently overridden
     unless the user explicitly specifies them by passing the arguments through
     \code{dots}. Such arguments include: \code{MemoryMatrix = FALSE}, 
     \code{print.level = 1}, \code{P9mix = 1}, \code{MemoryMatrix = FALSE}, 
     \code{print.level = 1}, \code{max.generations = 1000}, \code{project.path}, 
     \code{starting.values}, and --- in most cases --- 
     \code{boundary.enforcement = 1}. The arguments to \code{\link[rgenoud]{genoud}} 
     that remain at their defaults but you may want to seriously consider 
     tweaking are \code{pop.size} and \code{wait.generations}.}
}
\details{
     SEFA, EFA, and CFA models all estimate the same population model but impose
     different restrictions on the model. If \code{restrictions} is unspecified,
     \code{Factanal} will create an object that inherits from class
     \code{"restrictions"} based on the responses the user gives to the pop-up
     menus. The vignette provides a step-by-step guide to navigating the pop-up
     menus; execute \code{vignette("FAiR")} to read it.
     \code{Factanal} will then impose different restrictions on the model, 
     depending on the inherited class of this object.

     The CFA model is perhaps the most straightforward in the sense that the
     user specifies that certain coefficients are pegged to particular values, and 
     \code{Factanal} estimates the values of the free parameters. These 
     restrictions can be specified via the \code{fixed} argument or left 
     unspecified in which case \code{Factanal} will prompt you to specify
     the restrictions via a pop-up menu. \code{Factanal} relies on the theorem
     in Howe (1955) to overcome rotational indeterminancy. Namely, the factors are
     scaled to have unit variance  and at least \eqn{factors - 1} coefficients per 
     factor are pegged to zero such that a technical rank condition is satisfied.
     This mechanism for eliminating rotational indeterminancy is somewhat more
     limited than the options that are available in some other software packages
     for factor analysis but easily generalizes to SEFA models.

     A SEFA model differs from a CFA model in that the analyst specifies 
     how many coefficients per factor take the value of zero, and the algorithm 
     estimates the locations of these zeros along with the values of the 
     corresponding free coefficients. It is also possible to estimate a mixed 
     SEFA model where some coefficients are fixed to zero (or another number)
     \emph{a priori} and the locations of the remaining zeros are estimated.
     A SEFA model requires that the Howe (1955) theorem be satisfied and at least
     one additional restriction is imposed. SEFA models are new to the literature and 
     more information about them can be found in Goodrich (2008).

     A EFA model specifies an arbitrary set of restrictions that are minimally
     necessary to extract factors and then a transformation of the factors should
     be obtained using \code{\link{Rotate}}. By default, the fitting function
     used by \code{Factanal} to estimate a EFA model is the same as that used in
     \code{\link[stats]{factanal}}. However, there is an alternative choice that
     estimates a EFA model via a CFA algorithm with the upper triangle of the 
     coefficient matrix filled with zeros. The results, when \code{method = "MLE"}
     should be the same --- up to a transformation of the factors --- but in 
     practice can differ if there are optimization failures. The default algorithm
     is considered more reliable at this point, but the alternative algorithm
     must be used if \code{method = "YWLS"}. Which of these two algorithms is
     used depends on the class of \code{restrictions}, which in the usual case
     that it is unspecified will result in a pop-up menu asking which algorithm
     to use.

     It is not necessary to provide starting values for the parameters, since
     there are methods for that purpose. See \code{\link{S4GenericsFAiR}}. But a
     matrix of starting values can be passed to through the \code{dots} to
     \code{\link[rgenoud]{genoud}}. This matrix should have rows equal to the
     \code{pop.size} argument in \code{\link[rgenoud]{genoud}} and columns equal
     to the number of free free parameters in the model, which corresponds to
     the \code{nvars} argument in \code{\link[rgenoud]{genoud}}. The order of
     the parameters / columns proceeds from the \dQuote{top} of the model to
     the \dQuote{bottom} as follows. First come the cells that comprise the
     upper triangle of the factor intercorrelation matrix at level two 
     (if there is more than one second-order factor). Next come the 
     free cells of the coefficient matrix at level two in row-major order. 
     Then come the free cells that comprise the upper triange of the factor
     intercorrelation matrix at level one (if a second-order model is not
     estimated). Next come the free cells of the coefficient matrix at level
     one in row-major order. Finally come the diagonal cells of the uniqueness
     matrix. Note that a parameter is free unless it is fixed \emph{a priori},
     which is to say that in SEFA models coefficients are considered free
     even if there is a possibility that the algorithm will bind them to zero 
     at the optimum, rendering them \dQuote{not free} for the purpose of counting 
     degrees of freedom.
}
\value{
  An object of formal S4 class \code{"FA"}, or in the case of two-level models
  an object of formal S4 class \code{"FA.general"} or \code{"FA.2ndorder"}.
}
\note{
  The underlying genetic algorithm will print a variety of output as it progresses.
  On Windows, you have to move the scrollbar periodically to flush the output
  to the screen. The output will look something like this
  \tabular{llllll}{
         0 \tab 1.0   \tab 1.0 \tab \dots \tab 1.0 \tab double\cr
         1 \tab 1.0   \tab 1.0 \tab \dots \tab 1.0 \tab double\cr
         \dots \tab \dots \tab \dots \tab \dots \tab \dots \tab \dots\cr
         437 \tab 1.0   \tab 1.0 \tab \dots \tab 1.0 \tab double\cr
       }
  The integer on the very left indicates the generation number. If it appears to
  skip one or more generations, that signifies that the best individual in the
  \dQuote{missing} generation was no better than the best individual in the 
  previous generation. The sequence of ones indicates that various constraints
  are being satisfied by the best individual in the generation. Some of these
  constraints are hard-coded, some are added by the choices the user makes. The
  curious are referred to the source code, but for the most part users need not
  worry about them provided they are \eqn{1.0}. If any but the last are not \eqn{1.0}
  after the first few generations, there is a problem because no individual is 
  satisfying all the constraints. The last number is a double-precision number,
  typically the log-likelihood. This number will increase, sometimes painfully
  slowly, sometimes intermittently, over the generations since the log-likelihood
  is being maximized subject to the aforementioned constraints.
}
\section{Warning}{
  Yates' (1987 p.229) weighted least squares criterion has never received much 
  scrutiny. The criterion (but not Yates' algorithm) is included in FAiR so that
  it can be fully evaluated. However, it is not scale invariant, it does not
  lend itself to calculating standard errors or test statistics, and in limited 
  testing seems prone to finding a solution that is geared more toward minimizing 
  the weights than minimizing the squared residuals.
}
\references{
  Barthlomew, D. J.  and Knott, M. (1990) \emph{Latent Variable Analysis
    and Factor Analysis.} Second Edition, Arnold.

  Beauducel, A. (2007) In spite of indeterminancy, many common factor score
    estimates yield an identical reproduced covariance matrix.
    \emph{Psychometrika}, \bold{72}, 437--441.

  Goodrich, B. (2008) SEFAiR So Far. Unpublished manuscript linked at
    \url{http://wiki.r-project.org/rwiki/doku.php?id=packages:cran:fair#to_paper_s_about_the_ideas_in_fair}.

  Smith, G. A. and Stanley G. (1983)
  Clocking \eqn{g}: relating intelligence and measures of timed
  performance. \emph{Intelligence}, \bold{7}, 353--368.

  Venables, W. N. and Ripley, B. D. (2002)
  \emph{Modern Applied Statistics with S.} Fourth edition.  Springer.

  Yates, A. (1987) \emph{Multivariate Exploratory Data Analysis:
  A Perspective on Exploratory Factor Analysis.} 
  State University of New York Press.
}
\author{Ben Goodrich \url{http://wiki.r-project.org/rwiki/doku.php?id=packages:cran:fair}}
\seealso{ \code{\link{Rotate}} and \code{\link[stats]{factanal}} }
\examples{
## Example from Venables and Ripley (2002, p. 323)
## Previously from Bartholomew and Knott  (1999, p. 68--72)
## Originally from Smith and Stanley (1983)

data(ability.cov)
print(ability.cov)

if(TRUE){ # NOTE: One would usually not bother with this block. It just makes the
          # example go quickly and without user intervention on the pop-up menus.
starts1 <- c(0.4551693481819578, 
             0.5893203083906567, 
             0.2182044732474321, 
             0.7694294930481663,
             0.0526383747875095, 
             0.3334323600411430)
starts1 <- matrix(starts1, nrow = 1)

example1 <- new("restrictions.factanal", factors = 2L, nvars = 6L,
                Domains = cbind(sqrt(.Machine$double.eps), rep(1, 6)),
                model = "EFA", method = "MLE", dof = 4L, fast = FALSE)
}

# 'restrictions' and 'starting.values' would typically be left unspecified!
efa <- Factanal(covmat = ability.cov, factors = 2, model = "EFA",
                restrictions = example1, starting.values = starts1)
show(efa)
summary(efa)

# 'criteria' would typically be left unspecified!
efa.rotated <- Rotate(efa, criteria = list("phi"))
summary(efa.rotated)

if(TRUE){ # NOTE: One would usually not bother with this block. It just makes the
          # example go quickly and without user intervention on the pop-up menus.
starts2 <- c(4.46294498156615e-01,
             4.67036349420035e-01,
             6.42220238211291e-01,
             8.88564379236454e-01,
             4.77779639176941e-01,
            -7.13405536379741e-02,
            -9.47782525342137e-08,
             4.04993872375487e-01,
            -1.04604290549591e-08,
            -9.44950629176182e-03,
             2.63078925240678e-04,
             9.38038168787216e-01,
             8.43618801925473e-01,
             4.49024212016027e-01,
             5.87550265675745e-01,
             2.17850254355888e-01,
             7.71724777627142e-01,
             1.20084009542348e-01,
             2.88308011310065e-01)

starts2 <- matrix(starts2, nrow = 1)

Domains <- cbind(-1, 1)
Domains <- rbind(Domains, cbind(-1.5, rep(1.5, 12)))
Domains <- rbind(Domains, cbind(0, rep(1, 6)))
fixed   <- matrix(NA_real_, nrow = 6, ncol = 2)
fix_beta_args <- as.list(formals(FAiR:::FAiR_fix_coefficients))
fix_beta_args$zeros <- c(2,2)
beta_select <- c(FALSE, rep(TRUE, length(fixed)), rep(FALSE, nrow(fixed)))
beta_list <- list(beta = fixed, free = c(is.na(fixed)),
                  num_free = length(fixed), select = beta_select,
                  fix_beta_args = fix_beta_args)
Theta2_list <- list(Theta2 = diag(nrow(fixed)), 
                    select = c(rep(FALSE, length(fixed) + 1),
                               rep(TRUE, nrow(fixed))))
Phi <- diag(c(0.5, 0.5))
example2 <- new("restrictions.1storder", factors = c(2L, 0L),
                Domains = Domains, nvars = nrow(Domains), 
                model = "SEFA", method = "MLE", dof = 6L,
                Phi = Phi, beta = beta_list, Theta2 = Theta2_list,
                criteria = list(llik = FAiR:::FAiR_criterion_llik))
}

# 'restrictions' and 'starting.values' would typically be left unspecified!
sefa <- Factanal(covmat = ability.cov, factors = 2, model = "SEFA",
                 restrictions = example2, starting.values = starts2)
show(sefa)
summary(sefa)

stuff <- list() # output list for various methods, also works on efa and efa.rotated
stuff$model.matrix <- model.matrix(sefa) # sample correlation matrix
stuff$fitted <- fitted(sefa) # reproduced correlation with communalities on diagonal
stuff$residuals <- residuals(sefa) # difference between model.matrix and fitted
stuff$rstandard <- rstandard(sefa) # residual matrix rescaled to a correlation matrix
stuff$weights <- weights(sefa) # (scaled) approximate weights for residuals
stuff$influence <- influence(sefa) # weights * residuals
stuff$logLik <- logLik(sefa) # log-likelihood
stuff$BIC <- BIC(sefa) # BIC
stuff$profile <- profile(sefa) # profile plots of non-free parameters
plot(sefa)  # advanced Scree plot
pairs(sefa) # Thurstone-style plot
}
\keyword{ multivariate }
\keyword{ models }
