\name{mapping_rule}
\alias{mapping_rule}
\title{Default Mapping Rule}
\description{
  It is not generally necessary to understand this help page if one
  merely wants to estimate a factor analysis model and this function
  should never be called directly. However, it is used as the default
  function to squash some coefficients to zero during semi-exploratory 
  factor analysis and thus needs to be documented. Also, if one were 
  to specify a different mapping rule function in the call to 
  \code{\link{make_restrictions}}, then that function would need to 
  have some of the same arguments as this one.
}
\usage{
mapping_rule(coefs, cormat, zeros = rep(ncol(coefs), ncol(coefs)),
row_complexity = NA_integer_, quasi_Yates = FALSE, weak_Thurstone = FALSE, 
Butler = FALSE,  viral = FALSE, mve = FALSE, communality = FALSE)
}
\arguments{
  \item{coefs}{A primary pattern matrix with all cells filled in }
  \item{cormat}{A correlation matrix among primary factors with all 
    cells filled in}
  \item{zeros}{An integer vector whose length is the same as the
    number of factors that indicates the requisite number of zeros
    in the corresponding column of the primary pattern matrix}
  \item{row\_complexity}{An integer vector, either of length one
    or the same length as the number of rows in \code{coefs}, which
    indicates the requisite number of \emph{non-zero} coefficients
    in that row of the primary pattern matrix. If a scalar, this
    complexity is used for all rows of the primary pattern matrix.
    If \code{NA_integer_} (the default) this row-wise mapping rule 
    is not enforced }
  \item{quasi\_Yates}{A logical indicating whether to enforce a 
    mapping rule that is inspired by Yates (1987) to produce
    cohyperplanarity; see Details}
  \item{weak\_Thurstone}{ A logical indicating whether to enforce
    a mapping rule that is based on a minimal version of Thurstone's
    (1947) second and third rules of thumb for identifying the
    structure; see Details }
  \item{Butler}{ A logical indicating whether to enforce a mapping
    a mapping rule that is inspired by Butler (1969) to produce
    a unit complexity basis for the common factor space; see 
    Details}
  \item{viral}{ A logical indicating whether to enforce a mapping
    rule that puts some outcomes (namely double the number of
    factors) into two hyperplanes in common factor space; see
    Details }
  \item{mve}{ A logical that must be \code{FALSE} because this
    mapping rule is (perhaps permanently) disabled }
  \item{communality}{ A logical indicating whether to enforce a
    mapping rule that tries to put one zero per factor within
    a row of the primary pattern matrix that has a relatively
    large communality; see Details}
}
\details{
  The \dQuote{mapping rule} is a deterministic function that takes
  a matrix with an insufficient number of exact zeros (often no
  zeros at all) and transforms it rather minimally to a matrix with
  the requisite number of exact zeros in the specified columns and / 
  or rows. Hence, when the transformed matrix is returned, it should
  have enough zeros in the proper locations to potentially satisfy
  some theorem on model identification, typically 
  \enc{Reiersøl's}{Reiersol's} (1950). It is not possible to verify
  that all aspects of \enc{Reiersøl's}{Reiersol's} (1950) theorem
  are satisfied, but another internal function verifies as much as
  it can. During the (genetic) optimization, the discrepancy function 
  is called \emph{after} the mapping rule is enforced and after the
  primary pattern matrix is checked for compliance with
  \enc{Reiersøl's}{Reiersol's} (1950) theorem. Hence, optimization
  is intended to find the optimal estimates, subject to the 
  restriction that the model is identified (with the help of the
  mapping rule).

  By default, this function is copied into the \code{mapping_rule}
  slot in the call to \code{\link{make_restrictions}} and its
  \code{\link{formals}} are adjusted to correspond to the specified
  mapping rule. You may wish to do so manually in rare circumstances.
  If a different function is used to enforce a mapping rule, it
  should also have \code{coefs}, \code{cormat}, and \code{zeros}
  arguments as documented above.

  The vignette has more details on \enc{Reiersøl's}{Reiersol's} (1950) 
  theorem and defines these mapping rules in symbols. To describe
  them briefly, the default \dQuote{mapping rule} simply loops through the
  \code{zeros} argument and squashes the smallest \code{zeros[p]}
  cells (by magnitude) in the \eqn{p}th column of \code{coefs} to
  zero and returns the resulting matrix. This is the behavior
  depicted in the first example and can be brought about by leaving
  all the arguments that have default values at their defaults.

  Only one mapping rule should be requested. The execption is that
  the default mapping rule is often called after a non-default
  mapping rule has been applied to make sure that that the \eqn{p}th
  column has \code{zeros[p]} zeros in it.

  The \code{weak_Thurstone} mapping rule is similar but has an
  additional provision that no row of \code{coefs} may contain
  more than one (exact) zero. This mapping rule can
  be seen as minimally satisfying Thurstone's (1947) second and
  third rules for identification.

  The \code{row_complexity} mapping rule is simply applied to
  the rows of the primary pattern matrix rather than the columns.
  Note that the \dQuote{complexity} of an outcome is its number of 
  \emph{non-zero} coefficients. This mapping rule could loop over 
  \code{row_complexity} (it actually utilizes \code{\link{apply}}) 
  and squashes all but the largest (by magnitude) 
  \code{row_complexity[j]} coefficients in the \eqn{j}th row of the 
  reference structure matrix to zero and then rescales the result 
  back to the primary pattern matrix. Then the default mapping rule 
  is called, which may do nothing if there are now enough zeros in 
  each column of \code{coefs} and the result is returned.

  The \code{Butler} mapping rule concentrates the zeros within rows
  so that each factor correspons to an outcome of complexity one.
  Then, the default mapping rule is called to obtain more zeros
  for the \eqn{p}th factor.

  The \code{viral} mapping rule is less drastic but concentrates
  two zeros into each of \eqn{2 * factors} rows of \code{coefs}
  and is possibly useful if the number of factors is large. Then,
  the default mapping rule is called to obtain more zeros for the
  the \eqn{p}th factor.

  The other mapping rules are more advanced and rely on the factor
  contribution matrix which is the element-by-element product of
  the the pattern and structure matrices. The 
  \code{\link{rowSums}} of the factor contribution matrix is
  a vector of communalities. The \code{communality} mapping rule places
  one zero per factor in a row of \code{coefs} corresponding to
  a row of the factor contribution matrix with a high ratio of
  its arithmetic mean to its geometric mean. Then, the default
  mapping rule is called to obtain a sufficient number of zeros
  for the \eqn{p}th factor.

  The \code{quasi_Yates} mapping rule places zeros in rows of
  \code{coefs} corresponding to rows of the factor contribution
  matrix with large differences between columns. The idea is to
  place zeros in rows where one factor is weak and another factor
  is strong in terms of explaining the variance in the corresponding
  outcome variable. This mapping rule is intended to achieve 
  \dQuote{cohyperplanarity} as discussed in Yates (1987), albeit
  not in reference to SEFA.
}
\value{
  A matrix, specifically the primary pattern matrix with the
  requisite number of zeros in its columns and / or rows.
}
\author{ Ben Goodrich }
\seealso{\code{\link{make_restrictions}}, \code{\link{parameter.coef.SEFA-class}}}
\references{
  \enc{Reiersøl}{Reiersol}, O. (1950) On the Identifiability of Parameters in Thurstone's 
     Multiple Factor Analysis. \emph{Psychometrika}, \bold{15}, 121--149.
}
\examples{
  ## This is just a demo; you should NOT call this function directly

  cormat <- diag(2)        # factor intercorrelation matrix
  coefs  <- matrix(rnorm(20), nrow = 10, ncol = 2) # primary pattern matrix
  zeros <- c(2, 2)         # we require two zeros per factor

  any(coefs == 0)          # FALSE

  ## Default mapping rule
  coefs_default <- mapping_rule(coefs, cormat, zeros)
  colSums(coefs_default == 0)      # c(2, 2)
  # Now the 2 smallest coefficients in each column are squashed to zero
  print(cbind(coefs, NA, coefs_default), digits = 3)

  ## row_complexity mapping rule
  coefs_row <- mapping_rule(coefs, cormat, zeros, row_complexity = 1)
  colSums(coefs_default == 0)      # at least two per factor
  # Now the smaller coefficient in each row is squashed to zero
  print(cbind(coefs, NA, coefs_row), digits = 3)

  ## The other mapping rules are sort of useless in the two factor case

}
\keyword{ manip }
