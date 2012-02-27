\name{model_comparison}
\alias{model_comparison}
\alias{paired_comparison}
\title{ Compare Factor Analysis Models }
\description{
  These functions produce the usual model comparison statistics for factor analysis
  models.
}
\usage{
model_comparison(..., correction = c("swain", "bartlett", "none"), 
                 conf.level = .9, nsim = 1001)
paired_comparison(M_0, M_1)
}
\arguments{
  \item{...}{objects of \code{\link{FA-class}} produced by \code{\link{Factanal}}}
  \item{correction}{ character string indicating what correction to use}
  \item{conf.level}{ confidence interval for the RMSEA statistic}
  \item{nsim}{ number of simulations for the nonparametric tests, see Details}
  \item{M\_0}{object of \code{\link{FA-class}} produced by \code{\link{Factanal}} that nests
    \code{M_1}}
  \item{M\_1}{object of \code{\link{FA-class}} produced by \code{\link{Factanal}} that is
    nested within \code{M_0}}
}
\details{
For exactly two nested models, \code{paired_comparison} performs the simple version of the
test recommened in Satorra and Bentler (2000); however, it is up to the user to verify that
\code{M_1} is nested within \code{M_0}.

Any number of objects of \code{\link{FA-class}} that are produced by \code{\link{Factanal}}
can be passed to \code{model_comparison} and a wide variety of statistic tests and fit
indices will be calculated. The exact behavior heavily depends on how the model was estimated
and in the case of traditional maximum likelihood estimation also depends on the \code{correction} 
argument.

If \code{correction = "swain"} (the default), the maximum likelihood test statistic is scaled by one of the
correction factors in Swain (1975) that has been recommended in Herzog, Boomsma, and Reinecke (2007) 
and in Herzog and Boomsma (forthcoming) and is based on \url{http://www.ppsw.rug.nl/~boomsma/swain.R}
Users should refer to these works for details, simulation results, and in publications making use of 
this Swain correction. If \code{correction = "bartlett"}, the correction factor recommended in Bartlett 
(1950), which is only strictly appropriate for exploratory factor analysis and has been implemented in 
\code{\link{factanal}} for a long time. If \code{correction = "none"}, then no correction factor is 
utilized, which is also the behavior for models that do not use the traditional maximum likelihood 
discrepancy function. If the ADF discrepancy function is used (or one of its special cases), the primary test statistic is that advocated in Yuan and Bentler (1998) but the test in equation 2.20b of Browne (1984) is also calculated.

The (primary) test statistic is then used in the root mean squared error of approximation (RMSEA) (see
Steiger and Lind 1980) to conduct a test of \dQuote{close fit}, namely that the true RMSEA is less than 
\eqn{0.05}. Confidence intervals are also reported and depend on the value of \code{conf.level}. The RMSEA 
is in turn used to calculate Steiger's (1989) \eqn{\gamma}{gamma} index. In the maximum likelihood case,
both of these are affected by the \code{correction} factor.

If the traditional maximum likelihood discrepancy function is used, then the BIC and
SIC (Stochastic Information Criterion, see Preacher 2006 and Preacher, Cai, and MacCallum 2007) are
calculated. These information criteria can be used to compare nonnested models and in both cases,
smaller is better.

Finally, several model comparison statistics are calculated, largely based on the \code{summary.sem}
function in the \pkg{sem} package. Most of these statistics are discussed in Bollen (1989). These are
\tabular{ll}{
  \strong{List element} \tab \strong{Reference}\cr
  \code{McDonald} \tab McDonald's (1989) Centrality Index\cr
  \code{GFI} \tab \enc{Jöreskog's}{Joreskog's} and \enc{Sorböm's}{Sorbom's} (1981) Goodness of Fit Index\cr
  \code{AGFI} \tab \enc{Jöreskog's}{Joreskog's} and \enc{Sorböm's}{Sorbom's} (1981) Adjusted Goodness of Fit 
     Index\cr
  \code{SRMR} \tab Bentler's (1995) Standardized Root Mean-squared Residual\cr
  \code{TLI} \tab Tucker and Lewis (1973) Index\cr
  \code{CFI} \tab Bentler's (1995) Comparative Fit Index\cr
  \code{NFI} \tab Bentler and Bonett's (1980) Normalized Fit Index\cr
  \code{NNFI} \tab Bentler and Bonett's (1980) Nonnormalized Fit Index\cr
}
}
\section{Warning}{
If Yates' weighted least squares discrepancy function is used, the test statistic is not
strictly valid.
}
\value{
  \code{paired_comparison} produces an object of S3 class "htest"; \code{model_comparison}
  produces a list with the following elements:
  \item{restrictions}{the \code{restrictions} object for the model}
  \item{exact\_fit}{a list of one or more objects, usually of S3 class "htest", indicating
    the results of the associated test(s) of exact fit}
  \item{close\_fit}{a list of one or more objects, usually of S3 class "htest", indicating
    the results of the associated test or measure of \dQuote{close} fit}
  \item{fit\_indices}{a list of numeric fit indices}
  \item{infocriteria}{in the case of maximum likelihood estimation, a list of the 
    information criteria that were calculated}
}
\author{Ben Goodrich with major contributions by Anne Boomsma, John Fox, and Walter Herzog}
\seealso{ \code{\link{Factanal}} }
\references{
    Bentler, P.M. (1995), \emph{EQS structural equations program manual}. Encino, CA: Multivariate Software. 

    Bentler, P.M., & Bonett, D.G. (1980), \dQuote{Significance tests and goodness of fit in the analysis of 
      covariance structures}. \emph{Psychological Bulletin}, \bold{88}, 588--606.

    Browne, M.W. (1984), \dQuote{Asymptotically distribution-free methods for the analysis of covariance
      structures}, \emph{British Journal of Mathematical and Statistical Psychology}, \bold{37}, 62--83.

    Bollen, K. A. (1989) \emph{Structural Equations With Latent Variables}. Wiley.

    Herzog, W., and Boomsma, A. (forthcoming), \dQuote{Small-Sample Robust Estimators of
       Noncentrality-Based and Incremental Model Fit} \emph{Structural Equation Modeling}.

    Herzog, W., Boomsma, A., and Reinecke, S. (2007), \dQuote{The model-size effect on traditional 
      and modified tests of covariance structures}. \emph{Structural Equation Modeling}, \bold{14}, 
      361--390.

    Hotelling, H. (1931), \dQuote{The generalization of Student's ratio}, \emph{Annals of Mathematical 
      Statistics}, \bold{2}, 360--378.

    \enc{Jöreskog}{Joreskog}, K. G., and \enc{Sorböm}{Sorbom}, D. (1981). \emph{LISREL V: Analysis
     of linear structural relations by the method of maximum likelihood}. Chicago: Internationa
     Educational Services.

    McDonald, R.P. (1989), \dQuote{An index of goodness-of-fit based on noncentrality}, 
       \emph{Journal of Classification}, \bold{6}, 97--103.

    Preacher, K.J. (2006), \dQuote{Quantifying Parsimony in Structural Equation Modeling},
      \emph{Multivariate Behavioral Research} \bold{41}, 227--259.

    Preacher, K.J., Cai, L., and MacCallum, R.C. (2007), \dQuote{Alternatives to traditional
    model comparison strategies for covariance structure models.} in
    \emph{Modeling Contextual Effects in Longitudinal Studies}, eds.
    Little, T.D., Bovaird, J.A., and Card, N.A. Psychology Press.

    Satorra, A and Bentler, P.M. (2001), \dQuote{A scaled difference chi-square test statistic for 
      moment structure analysis,} \emph{Psychometrika}, \bold{66}, 507--514.

    Steiger, J.H. and Lind, J.C. (1980), \dQuote{Statistically based tests for the number of common factors}
      Paper presented at the annual meeting of the Psychometric Society, Iowa City, IA.

    Steiger, J.H. (1989), EzPATH: A supplementary module for SYSTAT and SYGRAPH. 
       Evanston, IL: SYSTAT.

    Swain, A.J. (1975). \emph{Analysis of parametric structures for variance matrices}. Unpublished doctoral
       dissertation, Department of Statistics, University of Adelaide, Australia.

    Tucker, L. R, and Lewis, C. (1973), \dQuote{A reliability coefficient for maximum likelihood factor 
      analysis}. \emph{Psychometrika}, \bold{38}, 1--10.
}
\examples{
## See example in Factanal()
}
\keyword{ multivariate }
\keyword{ models }
