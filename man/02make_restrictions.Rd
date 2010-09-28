\name{make_restrictions}
\docType{methods}
\alias{make_restrictions}
\alias{make_restrictions-methods}
\alias{make_restrictions,manifest.basic,ANY,ANY,ANY,ANY,ANY-method}
\alias{make_restrictions,manifest.basic,missing,missing,missing,missing,missing-method}
\alias{make_restrictions,manifest.basic,parameter.scale,missing,missing,missing,missing-method}
\alias{make_restrictions,manifest.basic,missing,parameter.coef,missing,missing,missing-method}
\alias{make_restrictions,manifest.basic,missing,parameter.coef,missing,parameter.coef,missing-method}
\alias{make_restrictions,manifest.basic,missing,parameter.coef,missing,parameter.coef,parameter.cormat-method}
\alias{make_restrictions,manifest.basic,missing,parameter.coef,parameter.cormat,missing,missing-method}
\alias{make_restrictions,manifest.basic,parameter.scale,parameter.coef,missing,parameter.coef,missing-method}
\alias{make_restrictions,manifest.basic,parameter.scale,parameter.coef,missing,parameter.coef,parameter.cormat-method}
\alias{make_restrictions,manifest.basic,parameter.scale,parameter.coef,parameter.cormat,missing,missing-method}
\title{ Make an object of class "restrictions" }
\description{
  This function is intended for users and sets up the restrictions to be imposed on
  the factor analysis model, which includes specifying whether an exploratory, confirmatory,
  or semi-exploratory model is to be estimated and what discrepancy function to use. It
  is a prerequisite for calling \code{\link{Factanal}} and assumes you have already called
  \code{\link{make_manifest}}.

  However, it is not necessary to understand the function arguments in great detail because 
  most of the functionality is implemented via pop-up menus. The vignette provides the 
  \dQuote{substantive} documentation for this function and provides screenshots of the pop-up 
  menus; execute \code{vignette("FAiR")} to view it.   It is possible (though not recommended 
  in normal usage) to avoid the pop-up menus entirely, in which case it is necessary to 
  thoroughly understand the documentation here and in the vignette as well as the 
  \code{\link{restrictions-class}} definition.

  In technical terms, this S4 generic function is a constructor for objects of 
  \code{\link{restrictions-class}}. The arguments in the signature of the S4 generic function are:\cr
  \code{manifest}, \code{Omega}, \code{beta}, \code{Phi}, \code{Delta}, and \code{Xi}\cr
  The first S4 method defined immediately below is the one indended for the vast majority
  of usage. It \emph{requires} \code{manifest} to be specified, \emph{forbids} these Greek letters
  from being specified, and has a few optional arguments that can be specified to modify the
  default behavior.
}
\usage{
\S4method{make_restrictions}{manifest.basic,missing,missing,missing,missing,missing}(
manifest, factors = NULL, model = c("SEFA", "EFA", "CFA"), 
discrepancy = "default", nl_1 = NULL, nl_2 = NULL)

# USE THE ABOVE METHOD IN MOST CASES!!! 
# USE ONE OF THE FOLLOWING METHODS IN SIMULATIONS, ETC.

\S4method{make_restrictions}{manifest.basic,parameter.scale,missing,missing,missing,missing}(
manifest, Omega, discrepancy = "default")

# Above method implies zero common factors

\S4method{make_restrictions}{manifest.basic,missing,parameter.coef,missing,missing,missing}(
manifest, beta, discrepancy = "default")

# Above method implies EFA with orthogonal factors

\S4method{make_restrictions}{manifest.basic,missing,parameter.coef,parameter.cormat,missing,missing}(
manifest, beta, Phi, discrepancy = "default", 
criteria = list(), methodArgs = list())

# Adjacent methods imply correlated factors

\S4method{make_restrictions}{manifest.basic,parameter.scale,parameter.coef,parameter.cormat,missing,missing}(
manifest, Omega, beta, Phi, discrepancy = "default", 
criteria = list(), methodArgs = list())

\S4method{make_restrictions}{manifest.basic,missing,parameter.coef,missing,parameter.coef,missing}(
manifest, beta, Delta, discrepancy = "default", 
criteria = list(), methodArgs = list())

# Adjacent methods imply one second-order factor

\S4method{make_restrictions}{manifest.basic,parameter.scale,parameter.coef,missing,parameter.coef,missing}(
manifest, Omega, beta, Delta, discrepancy = "default", 
criteria = list(), methodArgs = list())

\S4method{make_restrictions}{manifest.basic,missing,parameter.coef,missing,parameter.coef,parameter.cormat}(
manifest, beta, Delta, Xi, discrepancy = "default", 
criteria = list(), methodArgs = list())

# Adjacent methods imply multiple second-order factors

\S4method{make_restrictions}{manifest.basic,parameter.scale,parameter.coef,missing,parameter.coef,parameter.cormat}(
manifest, Omega, beta, Delta, Xi, discrepancy = "default", 
criteria = list(), methodArgs = list())
}
\arguments{
  \item{manifest}{a \emph{required} object that inherits from \code{\link{manifest.basic-class}}}
  \item{Omega}{usually missing or else an object of \code{\link{parameter.scale-class}}}
  \item{beta}{usually missing or else an object that inherits from \code{\link{parameter.coef-class}}}
  \item{Phi}{usually missing or else an object of \code{\link{parameter.cormat-class}}}
  \item{Delta}{usually missing or else an object that inherits from \code{\link{parameter.coef-class}}}
  \item{Xi}{usally missing or else an object of \code{\link{parameter.cormat-class}}}
  \item{factors}{either \code{NULL} or a vector of length one or two indicating the number of
     factors at level one and level two}
  \item{model}{a character string among \code{"SEFA"} (default), \code{"CFA"}, and \code{"EFA"}
     indicating whether to estimate a semi-exploratory, confirmatory, or exploratory model}
  \item{discrepancy}{a character string among \code{"default"}, \code{"MLE"}, \code{"ADF"}, 
       \code{"ELLIPTICAL"}, \code{"HK"}, \code{"SHK"}, and \code{"YWLS"} indicating which
       discrepancy function to use. The default behavior is \code{"ADF"} if possible, otherwise
       \code{"MLE"}. See the Note section for details.}
  \item{nl\_1}{either \code{NULL} or a function with an argument called \code{x} that imposes
     nonlinear exact restrictions on some cells of the primary pattern matrix at level one. See the 
     Details section.}
  \item{nl\_2}{either \code{NULL} or a function with an argument called \code{x} that imposes
     nonlinear exact restrictions on some cells of the primary pattern matrix at level two. See the 
     Details section.}
  \item{criteria}{either \code{NULL} or a list (possibly of length zero) of functions or
     character strings naming functions to be used as criteria in the lexical optimization process.
     See the Note section.}
  \item{methodArgs}{a list (possibly of length zero) of required arguments to the functions
     listed in the \code{criteria} argument. See the Note section.}
}
\details{
Let the factor analysis model in the population be
\deqn{\Sigma = \Omega(\beta\Phi\beta^\prime + \Theta)\Omega}{Sigma = Omega(beta Phi beta' + Theta)Omega}
where \eqn{\Sigma}{Sigma} is the covariance matrix among outcome variables, \eqn{\Omega}{Omega}
is a diagonal matrix of standard deviations of the manifest variables, \eqn{\beta}{beta} is the 
primary pattern matrix (calibrated to standardized variables) with one column per factor, 
\eqn{\Phi}{Phi} is the correlation matrix among the primary factors, and \eqn{\Theta}{Theta} 
is the diagonal matrix of uniquenesses, which is fully determined by \eqn{\beta}{beta}, 
\eqn{\Phi}{Phi}, and the requirement that the matrix within parentheses has ones down its diagonal. 
Hence, \eqn{\beta\Phi\beta^\prime + \Theta}{beta Phi beta' + Theta} is the model's correlation matrix 
among outcome variables as a function of the factors. This parameterization is often attributed to
Cudeck (1989, equation 21) or to Krane and McDonald (1975) and will essentially yield the
correct results even if the model is fit to a correlation matrix.

The \code{make_restrictions} methods set up the right-hand side of the factor analysis model,
including the restrictions that are placed on the parameters. \pkg{FAiR} differs fundamentally
from other factor analysis software in that you can place inequality restrictions on functions
of multiple parameters, although this mechanism is not permitted during the factor extraction
stage of an exploratory factor analysis. The classes that inherit from 
\code{\link{restrictions-class}} typically have slots that are Greek letters, which are
objects that inherit from \code{\link{parameter-class}}. But the \code{\link{restrictions-class}}
has additional slots that contain other information about the model.

Exploratory factor analysis (EFA) requires a minimal set of restrictions and permits no additional
restrictions in the factor extraction stage. EFA preliminarily assumes the factors are orthogonal,
i.e. \eqn{\Phi = I}{Phi = I}, and either assumes that 
\eqn{\beta^\prime\Theta^{-1}\beta}{beta' Theta^-1 beta} is diagonal or that \eqn{\beta}{beta}
has all zeros in its upper triangle. However, after a transformation of the factors has been
chosen (see \code{\link{Rotate}}), EFA takes no strong position on how many or which cells of
\eqn{\beta}{beta} are (near) zero, unless target rotation or the simplimax criterion is used.

Confirmatory factor analysis (CFA) allows the user to choose which cells of \eqn{\beta}{beta} are 
exactly zero according to substantive theory and also permits other kinds of restrictions. 
Semi-exploratory factor analysis (SEFA) allows the user to choose how many cells in each column of 
\eqn{\beta}{beta} are zero but does not require the user to specify which cells are zero. SEFA also 
permits other kinds of restrictions, including constraining specific cells in \eqn{\beta}{beta} to 
be zero, as in a CFA.

For CFA and SEFA, it is optionally possible to estimate a \dQuote{two-level} model where the correlation 
matrix among first-order factors is a function of one or more second-order factors. Hence, let 
the second-order factor analysis model in the population be
\deqn{\Phi = \Delta\Xi\Delta^\prime + \Gamma}{Phi = Delta Xi Delta' + Gamma}
where \eqn{\Phi}{Phi} is the correlation matrix among first-order primary factors, \eqn{\Delta}{Delta} 
is the second-order primary pattern matrix with one column per factor, \eqn{\Xi}{Xi} is the correlation 
matrix among the second-order primary factors, and \eqn{\Gamma}{Gamma} is the diagonal matrix of 
second-order uniquenesses, which is fully determined by \eqn{\Delta}{Delta}, \eqn{\Xi}{Xi}, and the 
requirement that \eqn{\Phi}{Phi} has ones down its diagonal. Hence, in a two-level model,
\eqn{\Phi}{Phi} is restricted to be an exact function of \eqn{\Delta}{Delta} and \eqn{\Xi}{Xi} and
the parameters for the two levels are estimated simultaneously.

\strong{If you are at all unsure about what to do, use the first method listed in the Usage section, where 
\code{"manifest"} is the only required argument}. This method covers all the functionality of the other
methods and will walk you through all the necessary steps using pop-up menus. However, if you would like to 
impose any nonlinear exact restrictions on the cells of \eqn{\beta}{beta} or \eqn{\Delta}{Delta}, then you 
need to define such functions (whose first argument should be called \code{"x"}) in the global environment 
and specify them as the \code{nl_1} and / or \code{nl_2}. See also \code{\link{parameter.coef.nl-class}} for 
details.

The other methods will ask fewer questions via pop-up menus, and perhaps none at all (see the Note section 
below). Hence, they place more faith in the user to specify the additional arguments correctly. See the files
in the FAiR/tests directory for many examples of setting up models the hard way without resorting to the
pop-up menu system. You can stop reading this help page now if you are content with the pop-up menus. 

If only \code{"Omega"} is specified, an object of \code{\link{restrictions.independent-class}} will emerge and
has no common factors. \eqn{\Omega}{Omega} is the only free parameter to estimate, which is only useful in 
constructing a null model to calculate some fit indices for another model (and all of this would normally 
be handled automatically by \code{\link{model_comparison}} anyway). If \code{"Omega"} is specified along with
other arguments, then this object of \code{\link{parameter.scale-class}} is merely passed along and
will become a slot of the resulting object that inherits from \code{\link{restrictions-class}}.

If \code{"beta"} is specified, but not \code{"Phi"} or \code{"Delta"}, an object of 
\code{\link{restrictions.orthonormal-class}} will emerge, which is appropriate for EFA where the upper 
triangle of \eqn{\beta}{beta} contains all zeros. If one prefers the maximum-likelihood discrepancy 
function, there is a faster EFA algorithm that makes the assumption that
\eqn{\beta^\prime\Theta^{-1}\beta}{beta' Theta^-1 beta} is diagonal, which can be brought about by 
specifying \code{model = "EFA"} and \code{discrepancy = "MLE"} in the call to \code{make_restrictions} 
using the first method in the Usage section above. Doing so will produce an object of 
\code{\link{restrictions.factanal-class}}, whose other methods will reproduce the behavior of 
\code{\link{factanal}}.

If both \code{"beta"} and \code{"Phi"} are specified, but not \code{"Delta"}, an object of
\code{\link{restrictions.1storder-class}} will emerge, which is appropriate for CFA or SEFA
when there are no second-order factors, implying that \eqn{\Phi}{Phi} has no structure, other than
that required of a correlation matrix. Whether the model is a SEFA depends on whether
\code{"beta"} inherits from \code{\link{parameter.coef.SEFA-class}}.

If both \code{"beta"} and \code{"Delta"} are specified, but not \code{"Xi"}, an object of
\code{\link{restrictions.general-class}} will emerge, which is appropriate for CFA or SEFA
with exactly one second-order factor. The off-diagonals of \eqn{\Phi}{Phi} are 
\eqn{\Delta\Delta^\prime}{Delta Delta'}. Whether the model is SEFA
(at level one) depends on whether \code{"beta"} inherits from 
\code{\link{parameter.coef.SEFA-class}} and \code{"Delta"} cannot inherit from
\code{\link{parameter.coef.SEFA-class}}.

If \code{"beta"}, \code{"Delta"}, and \code{"Xi"} are specified, an object of
\code{\link{restrictions.2ndorder-class}} will emerge, which is appropriate for CFA or SEFA
with multiple second-order factors. Whether this two-level model is SEFA depends on whether \code{"beta"} 
and / or \code{"Delta"} inherit from \code{\link{parameter.coef.SEFA-class}}.
}
\note{
In order to enforce inequality restrictions, it is necessary to use a genetic algorithm for
lexical optimization (see \code{\link[rgenoud]{genoud}}). Lexcial optimization requires
a sequence of piecewise functions, each of which evaluates whether the restriction is satisfied,
and an ultimate continuous function, which in this case is a discrepancy function (see
Browne (1984)).

In normal usage, the constraints would be specified in response to pop-up menus. In abnormal usage, 
the pop-up menus can be circumvented by specifying the constraints via the \code{criteria} argument,
whose elements are user-definied functions or character strings naming constraint functions, which 
are discussed in more detail below. The primary situation where the pop-up menus must be avoided is 
when you cannot rely on user-intervention, such as (some types of) Monte Carlo simulations, replication
scripts, examples in these documentation pages, the files in FAiR/tests directory of the source
code, etc. Conversely, there are templates to do virtually anything without resorting to the pop-up
menus in the FAiR/tests directory and anyone who is interested should look there. Note also that
when doing simulations, you should almost certainly be setting the \code{seeds} argument for various
functions to \code{NULL} at \emph{every} opportunity or your simulations will be invalid.

The discrepancy function is named by the \code{"discrepancy"} argument and is automatically appended 
to serve as the ultimate lexical criterion and should not be explicitly specified in \code{criteria} 
under any circumstances.The default behavior of the \code{"discrepancy"} argument depends on the exact
class of the the \code{manifest} argument. If \code{manifest} inherits from \code{"manifest.data"} or 
is of class \code{"manifest.basic.userW"}, the default discrepancy function is Browne's (1984) 
asymptotically distribution-free (\code{"ADF"}) discrepancy function (in certain cases the weight 
matrix for the ADF discrepancy function is diagonal, which is often refered to as 
\dQuote{diagonally weighted least squares}). Otherwise, this argument will default to 
\code{"MLE"}, which assumes the manifest variables are distributed multivariate normal. 

Both \code{"ELLIPTICAL"} and \code{"HK"} are discrepancy functions that put further
structure on the (inverse) weight matrix used by the \code{"ADF"} discrepancy 
function, where the former assumes there is a single kurtosis parameter and the 
latter assumes a different kurtosis for each manifest variable. The \code{"SHK"} 
discrepancy function differs from the \code{"HK"} discrepancy function in that 
the sample kurtosis estimates are shrunk towards the median kurtosis estimate using 
the same idea as in the \code{\link[corpcor]{cov.shrink}} shrinkage estimator of 
variance in the suggested \pkg{corpcor} package.

If \code{method = "YWLS"}, Yates' (1987 p.229) weighted least squares \dQuote{discrepancy} 
function is used, which does not satisfy Browne's (1984) definition of a discrepancy
function. This function has never received much scrutiny and is included in \pkg{FAiR} so that
it can be fully evaluated. However, it does not lend itself to calculating standard errors or 
test statistics, and in limited testing seems prone to finding a solution that is geared more 
toward minimizing the weights than minimizing the squared residuals. Note that it also differs 
from the way in which \dQuote{weighted least squares} is typically used in the factor analysis 
literature and in particular the LISREL program. \dQuote{Weighted least squares} in LISREL
essentially corresponds to \code{"ADF"} in \pkg{FAiR}.

User-defined functions can be passed as elements of \code{criteria} and will be called in
the \code{\link{environment}} of the important but not exported function 
\code{FAiR:::FAiR_lexical_driver}, which is defined in FAiR/R/Utils.R. If \code{criteria} 
is a list that includes character strings, they should be names of one or more of the  
\dQuote{canned} functions in the following table. Some of these functions have additional 
arguments that can be specified as elements of \code{methodArgs} to avoid seeing a pop-up
menu asking for them. See the vignette for formal definitions of these restrictions and 
\code{methodArgs}.

\tabular{lll}{
  \strong{name} \tab \strong{methodArgs} \tab \strong{reminder of what function does}\cr
  \code{"ranks_rows_XXX"}  \tab \code{row_ranks} \tab row-wise ordering constraints\cr
  \code{"ranks_cols_XXX"}  \tab \code{col_ranks} \tab column-wise ordering constraints\cr
  \code{"indicators_XXX"} \tab \code{indicators} \tab designate which is the best indicator of a factor\cr
  \code{"evRF_XXX"}  \tab none \tab restrict effective variance of reference factors\cr
  \code{"evPF_XXX"}  \tab none \tab restrict effective variance of primary factors\cr  
  \code{"h2_over_FC_XXX"}  \tab none \tab communalities \eqn{>=} factor contributions\cr
  \code{"no_neg_suppressors_XXX"}  \tab \code{FC_threshold} \tab no negative suppressors\cr
  \code{"gv_XXX"}  \tab none \tab generalized variance of primary <= reference factors\cr
  \code{"distinguishability_XXX"}  \tab none \tab best indicators have no negative suppressors\cr
  \code{"cohyperplanarity_XXX"}  \tab none \tab hard to explain in this table\cr
  \code{"dist_cols_XXX"}  \tab \code{cutpoint} \tab minimum binary distance between columns\cr
  \code{"volume_1st"} \tab none \tab \code{det(R) >= det(S)}\cr
  \code{"block_XXX"}  \tab \code{blockers} \tab force coeffients to be non-zero in SEFA
}

   \code{"XXX"} can be either \code{"1st"} or \code{"2nd"} to indicate that the restriction
   should be applied at the first or second level of a model if there are two levels. However,
   \code{volume_1st} is defined only at the first level of the model. If, for example, one wanted to place 
   row-wise restrictions on the order of the factor contributions at level one and level two of a model, 
   then append \code{"_1st"} and \code{"_2nd"} to the corresponding elements of \code{methodArgs}.
}
\value{
  Returns an object that inherits from \code{\link{restrictions-class}}. This object would 
  then be passed to the \code{restrictions} argument of \code{\link{Factanal}}.
}
\author{Ben Goodrich}
\seealso{ \code{\link{Factanal}} and \code{\link{restrictions-class}} }
\references{
    Browne, M.W. (1984), Asymptotically distribution-free methods for the analysis of covariance
      structures, \emph{British Journal of Mathematical and Statistical Psychology}, \bold{37}, 62--83.

    Cudeck, R. (1989),  Analysis of correlation matrices using covariance structure models,
    \emph{Psychological Bulletin}, \bold{105} 317--327.

    Krane, W. R., and McDonald, R. P. (1978). Scale invariance and the factor analysis of covariance matrices. 
    \emph{British Journal of Mathematical and Statistical Psychology}, \bold{31}, 218–-228.

    Yates, A. (1987) \emph{Multivariate Exploratory Data Analysis:
    A Perspective on Exploratory Factor Analysis.} 
    State University of New York Press.
}
\examples{
man <- make_manifest(covmat = Harman74.cor)

\dontrun{
## Here is the easy way to set up a two-level SEFA model using pop-up menus
res <- make_restrictions(manifest = man, factors = c(5,2), model = "SEFA")
}

## Here is the hard way to set up a two-level SEFA model, 
## which eschews the pop-up menus; do NOT do this without good reason.
## There is an EFA example in ?Rotate
## There is a  CFA example in ?restrictions2RAM
## There is an example with equality restrictions in ?equality_restriction-class

factors <- c(5,2)
beta <- matrix(NA_real_, nrow = nrow(cormat(man)), ncol = factors[1])
rownames(beta) <- rownames(cormat(man))
free <- is.na(beta)
beta <- new("parameter.coef.SEFA", x = beta, free = free, num_free = sum(free))

Delta <- matrix(NA_real_, nrow = factors[1], ncol = factors[2])
rownames(Delta) <- paste("F",  1:factors[1], sep = "")
free  <- is.na(Delta)
Delta <- new("parameter.coef.SEFA", x = Delta, free = free, num_free = sum(free))

Xi <- diag(2)
free <- lower.tri(Xi)
## For fun, require the second-order factors to be positively correlated
uppers <- lowers <- Xi
lowers[2,1] <- 0
uppers[2,1] <- 1
Domains <- array(cbind(lowers, uppers), c(dim(Xi), 2))
Xi <- new("parameter.cormat", x = Xi, free = free, num_free = sum(free),
           Domains = Domains)

res <- make_restrictions(manifest = man, beta = beta, Delta = Delta, Xi = Xi)
show(res)
}
\keyword{ multivariate }
\keyword{ models }
\keyword{ methods }
