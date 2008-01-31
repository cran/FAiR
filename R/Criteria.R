#     This file is part of FAiR, a program to conduct Factor Analysis in R
#     Copyright 2008 Benjamin King Goodrich
# 
#     FAiR is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
# 
#     FAiR is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
# 
#     You should have received a copy of the GNU General Public License
#     along with FAiR.  If not, see <http://www.gnu.org/licenses/>.

## This file contains (some of) the functions that are used as lexical criteria

## NOTE: This file is meant to be read with 90 columns and 8 space tabs

## Ultimate criterions for Factanal()
FAiR_criterion_llik <- # log-likelihood
function() {
	# Both S and C are available in the calling environment
	# S is the sample correlation matrix, C is the estimate of Sigma
	C_inverse <- chol2inv(chol(C))
	llik <- -determinant(C)$modulus - crossprod(c(S), c(C_inverse))
	return(c(loglik = llik))
}

FAiR_criterion_YWLS <- # Yates' (1987) weighted least squares criterion
function() {
	# beta, S, and C are available in the calling environment
        # beta is the primary pattern matrix, C is the estimate of Sigma, 
	# S is the sample correlation matrix
	communalities <- rowSums(beta^2)
	diag(C) <- communalities
	w <- 1 - sweep(sweep(C^2, 1, communalities, FUN = "/", check.margin = FALSE),
                                  2, communalities, FUN = "/", check.margin = FALSE)
	return(c(NegWLS = -(sum(w * (S - C)^2))))
}

## Restriction criteria used exclusively in the factor extraction phase of SEFA and CFA
 # for the first four criteria, Phi is the correlation matrix among primary factors at 
 # level 1, Xi is the correlation matrix among primary factors at level 2 and C is the
 # model reproduced correlation matrix among outcomes with ones on the diagonal
 # all of these matrices are available in the calling environment
FAiR_criterion_evRF_2nd <- function() {
	ev_Xi <- det(Xi)^(1/ncol(Xi))
	ev_Phi <- det(Phi)^(1/ncol(Phi))
	ev_diff <- ev_Xi - ev_Phi
	if(ev_diff >= 0) return(1.0)
	else return(ev_diff)
}

FAiR_criterion_evPF_2nd <- function() {
	ev_RF <- det(cov2or(chol2inv(chol(Xi))))^(1/ncol(Xi))
	ev_Phi <- det(Phi)^(1/ncol(Phi))
	ev_diff <- ev_RF - ev_Phi
	if(ev_diff >= 0) return(1.0)
	else return(ev_diff)
}

FAiR_criterion_evRF_1st <- function() {
	ev_Psi <- det(cov2cor(chol2inv(chol(Phi))))^(1/ncol(Phi))
	ev_C <- det(C)^(1/ncol(C))
	ev_diff <- ev_Psi - ev_C
	if(ev_diff >= 0) return(1.0)
	else return(ev_diff)
}

FAiR_criterion_evPF_1st <- function() {
	ev_Phi <- det(Phi)^(1/ncol(Phi))
	ev_C <- det(C)^(1/ncol(C))
	ev_diff <- ev_Phi - ev_C
	if(ev_diff >= 0) return(1.0)
	else return(ev_diff)
}

FAiR_criterion_no_suppressors_1st <- # forces factor contributions to exceed a threshold
function(threshold = -0.01) {
	# Both beta, and Phi are available in the calling environment
	# beta is the coefficient matrix, Phi is the correlation matrix among factors
	# threshold is the minimum acceptable factor contribution
	Pi <- beta %*% Phi
	factor.contributions <- beta * Pi
	return(c(NS1 = mean(factor.contributions >= threshold)))
}

FAiR_criterion_no_suppressors_2nd <- # forces factor contributions to exceed a threshold
function(threshold = -0.01) {
	# Both Delta, and Xi are available in the calling environment
	# Delta is the coefficient matrix, Xi is the correlation matrix among factors
	# FC_threshold is the minimum acceptable factor contribution
	Gamma <- Delta %*% Xi
	factor.contributions <- Delta * Gamma
	return(c(NS1 = mean(factor.contributions >= threshold)))
}

FAiR_criterion_dets_2nd <- # forces primary factors to have less generalized variance
function() {               # than the corresponding reference factors
	# Xi is available in the calling environment and is the correlation matrix
	# among second-order primary factors
	D2 <- rowSums(backsolve(chol(Xi), diag(ncol(Xi)))^2)
	det_Xi <- det(Xi)
	out <-  det_Xi^2 - prod(D2) 
	if(out <= 0) return(c(det2 = 1.0))
	else         return(c(det2 = out))
}

FAiR_criterion_dets_1st <- # forces primary factors to have less generalized variance
function() {               # than the corresponding reference factors
	# Phi is available in the calling environment and is the correlation matrix
	# among first-order primary factors
	D2 <- rowSums(backsolve(chol(Phi), diag(ncol(Phi)))^2)
	det_Phi <- det(Phi)
	out <-  det_Phi^2 - prod(D2)
	if(out <= 0) return(c(det1 = 1.0))
	else         return(c(det1 = out))
}

FAiR_criterion_cohyperplanarity <- # forces tests within hyperplanes to have more 
function() {                       # effective variance than the battery as a whole
	# beta and C are available in the calling environment and are the primary pattern
	# matrix and the estimated correlation matrix (ones on the diagonal) respecitvely
	EV_C <- det(C)^(1/ncol(C))
	out  <- rep(NA, ncol(beta))
	for(i in 1:ncol(beta)) {
		mark <- beta[,i] == 0
		out[i] <- (det(C[mark,mark])^(1/sum(mark)) >= EV_C)
	}
	return(mean(out))
}

## Criteria used in Rotate()
## For all that follow, Upsilon is the reference structure matrix and Phi is the 
## correlation matrix among primary factors. Both are available in the calling environment

## Ultimate criterions
FAiR_criterion_phi <- # (log of) Thurstone's criterion, possibly with his scale factor
function(c = 1) {
	return(c(logphi = log(sum(exp(rowSums(log(Upsilon^2)/c))))))
}

FAiR_criterion_varphi <- # BG's generalization of Thurstone's criterion
function(weights = NULL) {
	# Each row of the reference structure matrix is sorted by magnitude
	sorted <- t(apply(Upsilon^2, 1, sort))

	# Start with a plain Thurstone's criterion
	varphi <- sum(exp(rowSums(log(sorted))))

	# Exclude a column, recalculate Thurstone's criterion, weight it, and cumulate
	if(is.null(weights)) for(i in 1:(ncol(sorted) - 1)) {
		weight <- max(0, 1 - max(sorted[,i])) # dynamic weights
		varphi <- varphi + weight    * sum(exp(rowSums(
						log(sorted[,-c(1:i),drop = FALSE]))))
	}
	else for(i in 1:(ncol(sorted) - 1)) {
		varphi <- varphi + weights[i] * sum(exp(rowSums(
						log(sorted[,-c(1:i),drop = FALSE]))))
	}
	return(c(logvarphi = log(varphi)))
}

FAiR_criterion_minimaximin <- # Another way of operationalizing simple structure
function() {	
	return(c(logmmm = log(max(apply(Upsilon^2, 1, min)))))
}

FAiR_criterion_LS <- # Loading Simplicity Index by Lorenzo-Seva (2004) p.50-1
function(eps = .Machine$double.eps, scale = 10,
           e = (1 / ncol(Upsilon) + eps)^(scale / ncol(Upsilon))) {
	# scale invariance implies it can be calculated on the reference structure matrix
	C <- diag(colSums(Upsilon^2))
	H <- diag(diag(sweep(Upsilon, 2, 1/diag(C), "*", FALSE) %*% t(Upsilon)))
	B <- sweep(sweep(Upsilon, 1, diag(H)^(-.5), "*", FALSE), 
				  2, diag(C)^(-.5), "*", FALSE)

	w <- mean( (B^2 + eps)^(scale * B^2) )

	LS <- (w - e) / (1 - e)
	if(is.nan(LS)) LS <- 1  # absolutely perfect clustering -> NaN
	return(-LS)
}

## Restriction criteria
FAiR_criterion_no_factor_collapse <- # forces effective variance of the primary factor
function(Phi, threshold = 0.25) {    # intercorrelation matrix to exceed a threshold
	# An EV of zero implies factor collapse so make the threshold 0.25 or something
	effective_variance <- det(Phi)^(1/ncol(Phi))
	if(effective_variance >= threshold) return(c(nocollapse = -1.0))
	else return(c(collapse = -effective_variance))
}

FAiR_criterion_limit_correlations <- # curtail pairwise correlations between primary
function(lower, upper) {             # factors
	# Can be used to make inter-factor correlations positive (i.e. lower = 0) or to
	# prevent factor fission (upper << 1). Note that no_factor_collapse and possibly
	# no_suppressors can be used to accomplish a similar goal much of the time.
	cors <- Phi[lower.tri(Phi)]
	toolow  <- cors < lower
	toohigh <- cors > upper
	inadmissable <- toolow | toohigh
	if(any(inadmissable)) {
		if(any(toolow))  x <- (cors[toolow]  - lower)^2
		else x <- 0

		if(any(toohigh)) y <- (cors[toohigh] - upper)^2
		else y <- 0

		return(c(outofbounds = sqrt(max(c(x,y)))))
	}
	else return(c(inbounds = -1.0))
}

FAiR_criterion_evRF <- function(threshold = NA) {
	ev_Psi <- det(cov2cor(chol2inv(chol(Phi))))^(1/ncol(Phi))
	ev_diff <- threshold - ev_Phi
	if(ev_diff <= 0) return(-1.0)
	else return(ev_diff)
}

FAiR_criterion_evPF <- function(threshold = NA) {
	ev_Phi <- det(Phi)^(1/ncol(Phi))
	ev_diff <- threshold - ev_Phi
	if(ev_diff <= 0) return(-1.0)
	else return(ev_diff)
}

FAiR_criterion_no_suppressors <- # forces factor contributions to exceed a threshold
function(threshold = -0.01) {
	RP <- A %*% rotmat_reference %*% chol2inv(chol(crossprod(rotmat_reference)))
	factor.contributions <- RP * Upsilon
	admissable <- mean(factor.contributions >= threshold)
	if(admissable < 1) return(c(suppressors = -min(factor.contributions)))
	else               return(c(NS1 = -1.0))
}

FAiR_criterion_dets <- # forces primary factors to have less generalized variance
function() {           # than the corresponding reference factors
	# Phi is available in the calling environment and is the correlation matrix
	# among first-order primary factors
	D2 <- rowSums(backsolve(chol(Phi), diag(ncol(Phi)))^2)
	det_Phi <- det(Phi)
	out <-  det_Phi^2 - prod(D2)
	if(out <= 0) return(c(det1 = -1.0))
	else         return(c(det1 =  out))
}

FAiR_criterion_positive_manifold <- # forces reference structure correlations to exceed
function(threshold = -0.1) {        # a threshold
	pm <- mean(Upsilon >= threshold)
	if(pm < 1) return(c(noPM = -min(Upsilon)))
	else       return(c(PM = -1.0))
}
