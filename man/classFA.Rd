\name{classFA}
\docType{class}
\alias{classFA}
\alias{FA-class}
\alias{FA.general-class}
\alias{FA.2ndorder-class}
\alias{BIC,FA-method}
\alias{confint,FA-method}
\alias{coef,FA-method}
\alias{logLik,FA-method}
\alias{plot,FA,ANY-method}
\alias{profile,FA-method}
\alias{show,FA-method}
\alias{summary,FA-method}
\alias{vcov,FA-method}
\alias{FA.general-class}
\alias{coef,FA.general-method}
\alias{FA.2ndorder-class}
\alias{coef,FA.2ndorder-method}
\title{Class "FA" and Classes that Inherit from "FA"}
\description{These classes encapsulate estimates from factor analysis 
  models. The \code{"FA.general"} and \code{"FA.2ndorder"} classes extend 
  class \code{"FA"} and have additional slots to contain estimates from a 
  second-order model with one general second-order factor and multiple 
  second-order factors respectively}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("FA", ...)}. However,
this use of \code{new("FA", ...)} is not recommended. Both 
\code{\link{Factanal}} and \code{\link{Rotate}} call \code{new("FA", ...)}
internally to create the appropriate object for analysis with the help
of the formal methods defined for \code{create_FAobject}.
}
\section{Slots}{
	 \describe{
    \item{\code{loadings}:}{Object of class \code{"array"}. This array
      has as many rows as there are outcome variables, as many columns
      as there are factors, and five shelves. Each shelf is a matrix
      and contains the primary pattern ('PP'), primary structure ('PS'), 
      reference pattern ('RP'), reference structure ('RS'), and factor 
      contribution ('FC') matrices respectively, which can be extracted via 
      \code{coef()} using the appropriate two-letter string for the 
      \code{matrix} argument that defaults to 'PP'.}
    \item{\code{correlations}:}{Object of class \code{"array"}. This array
      has as many rows and columns as there are factors and three shelves.
      Each shelf is a correlation matrix and contains the correlations
      among the primary factors ('PF'), among reference factors ('RF'), 
      and between each primary factor and the corresponding reference 
      factor ('PR'), which is a diagonal matrix. There is currently no
      special accessor function for these three matrices.}
    \item{\code{trans_mats}:}{Object of class \code{"array"}. This array
      has as many rows and columns as there are factors and two shelves.
      Both shelves are the identity matrix if semi-exploratory or 
      confirmatory factor analysis was used in the call to 
      \code{\link{Factanal}} and are the identity matrix prior to the
      call to \code{\link{Rotate}} in exploratory factor analysis. 
      \code{\link{Rotate}} fills these two shelves with the transformation
      matrix for the primary factors ('primary') and the  transformation
      matrix for the reference factors ('reference'). Both matrices have
      unit-length columns and the rotated primary pattern matrix is equal
      to the unrotated primary pattern matrix times the inverse of the
      transpose of the transformation matrix for the primary factors.}
    \item{\code{uniquenesses}:}{Object of class \code{"numeric"}. This
      numeric vector has as many elements as there are outcome variables
      and contains the unique variances (also called specific variances)
      of the factor analysis model.}
    \item{\code{restrictions}:}{Object of class \code{"restrictions"}. This
      slot contains information on the restrictions that were placed on 
      the factor analysis model while it was estimated.}
    \item{\code{vcov}:}{Object of class \code{"matrix"}. This matrix is
      a variance-covariance matrix and thus has as many rows and columns
      as there are parameters to be estimated. Note that it is calibrated
      to the sample covariance matrix, rather than the sample correlation
      matrix.}
    \item{\code{zstats}:}{Object of class \code{"list"}. A list of matrices
      containing the z-statistic for each estimated parameter under a null
      hypothesis that the parameter is zero. }
    \item{\code{scores}:}{Object of class \code{"matrix"}. If, in the 
      call to \code{\link{Factanal}} the user specifies that factor
      scores should be calculated, this matrix contains the factor scores
      and has as many rows as there are complete observations and as
      many columns as there are factors. If factor scores were not requested,
      then this matrix has one row, one column, and a NA\_real\_ in its 
      only cell.}
    \item{\code{manifest}:}{Object of class \code{"list"}. This list 
      contains information on the left-hand side of the factor analysis
      model. It has an element for the sample covariance matrix, the 
      sample correlation matrix, the number of observations in the sample,
      and possibly a matrix of standardized scores on the outcome variables. }
    \item{\code{rgenoud}:}{Object of class \code{"list"}. This is a list of 
      lists each containing what is returned by \code{\link[rgenoud]{genoud}}, 
      which is called internally by \code{\link{Factanal}} and 
      \code{\link{Rotate}}.}
    \item{\code{model}:}{Object of class \code{"character"}. This 
      character string indicates whether a SEFA, EFA, or CFA model was
      estimated.}
    \item{\code{method}:}{Object of class \code{"character"}. This 
      character string indicates whether maximum-likelihood or weighted
      least squares was used to estimate the factor analysis model. }
    \item{\code{call}:}{Object of class \code{"language"}.  This slot
      contains the call to \code{\link{Factanal}}.}
    \item{\code{seeds}:}{Object of class \code{"matrix"}. This matrix
      has two columns and either one or two rows. The first row contains
      the \code{unif.seed} used by \code{\link[rgenoud]{genoud}} and the
      \code{int.seed} in the call to \code{\link{Factanal}}. If 
      \code{\link{Rotate}} is used to transform the factors after they
      have been extracted via exploratory factor analysis, this matrix
      has a second row containing the \code{unif.seed} and 
      \code{int.seed} used in the call to \code{\link{Rotate}}.}
  }
}
\section{Methods}{
  \describe{
    \item{BIC}{\code{signature(object = "FA")}: Bayesian Information Criterion }
    \item{coef}{\code{signature(object = "FA")}: Extract coefficients.
      This method has two additional arguments, \code{matrix} and \code{level}.
      \code{matrix} should be one of \code{"PP"} (default) for the primary
      pattern matrix, \code{"PS"} for the primary structure matrix, 
      \code{"RP"} for the reference pattern matrix, \code{"RS"} for the
      reference structure matrix, or \code{"FC"} for the factor contribution
      matrix. \code{level} should be \code{1} (default) to extract 
      from level one of the model or \code{2} to extract from level two of
      the model (if the model has two levels).}
    \item{confint}{\code{signature(object = "FA")}: Extract confidence intervals.
      Does not work yet.}
    \item{logLik}{\code{signature(object = "FA")}: Extract log-likelihood,
      requires that the model be estimated by maximum likelihood.}
    \item{plot}{\code{signature(x = "FA", y = "ANY")}: Plots from \pkg{nFactors} }
    \item{profile}{\code{signature(fitted = "FA")}: Plot profile likelihoods
      for all fixed and zero coefficients in the primary pattern matrix. This
      method has two additional arguments, \code{param}, which is currently
      ignored and \code{delta} which indicates how far on either side of the
      constrained optimum to plot.}
    \item{show}{\code{signature(object = "FA")}: This method displays some
      information about the restrictions that were imposed on the model and
      prints a variety of model fit statistics (see note below). }
    \item{summary}{\code{signature(object = "FA")}: Gathers the coefficients,
      intercorrelations, and z-statistics, typically for showing  }
    \item{vcov}{\code{signature(object = "FA")}: Extracts the 
      variance-covariance matrix }
	 }
}
\note{The model comparison statistics printed by the \code{show} method are
largely resused from the \code{summary.sem} function in \pkg{sem}. Most of
these statistics are defined in Bollen (1989). The exception is the stochastic
information complexity (SIC), which is defined in Preacher, Cai, and MacCallum
(2007).}
\references{
    Bollen, K. A. (1989)
    \emph{Structural Equations With Latent Variables.} Wiley.

    Preacher, K.J., Cai, L., and MacCallum, R.C. ``Alternatives to traditional
    model comparison strategies for covariance structure models.'' in
    \emph{Modeling Contextual Effects in Longitudinal Studies}, eds.
    Little, T.D., Bovaird, J.A., and Card, N.A. Psychology Press.
}
\author{Ben Goodrich \url{http://wiki.r-project.org/rwiki/doku.php?id=packages:cran:fair}}
\seealso{
  \code{\link[stats4]{BIC}}, \code{\link[stats4]{coef}}, 
  \code{\link[stats4]{profile}}, \code{\link[stats4]{confint}}
  \code{\link[stats4]{logLik}}, \code{\link[stats4]{plot}}, \code{\link{show}},
  \code{\link{summary}}, \code{\link[stats4]{vcov}}, and \link{S3methodsFA}.
}
\examples{
showClass("FA")
}
\keyword{classes}
