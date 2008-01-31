#     This file is part of FAiR, a program to conduct Factor Analysis in R
#     Copyright 2008 Benjamin King Goodrich
#     Some portions of this code are Copyright 1995-2007 R Core Development Team,
#     Copyright 2007 Jens Oehlschl�gel, and Copyright John Fox, are licensed under the 
#     GPL version 2 or later, and are indicated so in the respective functions
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


## This file contains utility functions that are not meant to be accessed by users

## NOTE: This file is meant to be read with 90 columns

FAiR_parse <- # this function parses the inputs and returns a list
function(x, data = NULL, covmat = NULL, n.obs = NA, 
	subset, na.action, robust = FALSE, seeds) {
	## This function is slightly modified from stats:::factanal, which is 
        ## Copyright 1995-2007 R Core Development Team and licensed under GPL V2+
	set.seed(seeds)
	seeds <- get(".Random.seed", envir = .GlobalEnv,inherits = FALSE)
	na.act <- NULL
	have.x <- FALSE
	if (is.list(covmat)) {
		if (any(is.na(match(c("cov", "n.obs"), names(covmat))))) {
			stop("'covmat' is not a valid covariance list")
		}
		if(robust) stop("robust estimation of the covariance matrix is only",
				"possible if 'x' is specified rather than 'covmat'")
		if(!isTRUE(all.equal(covmat$cov, t(covmat$cov)))) {
			stop("covariance matrix is not symmetric")
		}
		if(!missing(x)) {
			stop("either pass 'x' (preferable) or 'covmat' but not both")
		}
		covmat$cor <- cov2cor(covmat$cov)
	}
	else if (is.matrix(covmat)) {
		if(nrow(covmat) != ncol(covmat)) stop("'covmat' must be symmetric")
		if(robust) stop("robust estimation of the covariance matrix is only",
				"possible if 'x' is specified rather than 'covmat'")
		if(!isTRUE(all.equal(covmat, t(covmat)))) {
			stop("covariance matrix is not symmetric")
		}
		if(!missing(x)) {
			stop("either pass 'x' (preferable) or 'covmat' but not both")
		}
		covmat <- list(cov = covmat, cor = cov2cor(covmat), n.obs = n.obs)
	}
	else if(!is.null(covmat)) {
		stop("covmat must be NULL, a list, or a covariance matrix")
	}
	else { # NULL covmat
        	if (missing(x)) stop("either 'x' or 'covmat' must be supplied")
		else have.x <- TRUE

        	if (inherits(x, "formula")) {
			mt <- terms(x, data = data)
			if (attr(mt, "response") > 0) {
				stop("response not allowed in formula")
			}
			attr(mt, "intercept") <- 0
			mf <- match.call(expand.dots = FALSE)
			names(mf)[names(mf) == "x"] <- "formula"
			mf[[1]] <- as.name("model.frame")
			mf[!names(mf) %in% c("", "formula", "data")] <- NULL
			mf <- eval.parent(mf)
			na.act <- attr(mf, "na.action")
			z <- model.matrix(mt, mf)
			if(all(sapply(mf, is.numeric))) {
				if(robust) {
					if(!require(robustbase)) { 
						stopifnot(require(MASS))
						covmat <- cov.mcd(z, cor = TRUE, nsamp = 
								"best", seed = seeds)
					}
					else {
						covmat <- covMcd(z, cor = TRUE, nsamp =
								"best", seed = seeds)
					}
				}
				else covmat <- cov.wt(z, cor = TRUE)
			}
			else {
				stop("ordinal variables not supported yet")
#				stopifnot(require(polycor))
				print("constructing correlation matrix among",
					"heterogenous variables which could take",
					"a while.")
				cormat <- hetcor(z, ML = TRUE, std.err = FALSE)
				covmat <- list(cov = NA_real_, n.obs = nrow(z), 
						cor = cormat)
			}
        	}
        	else {
			z <- as.matrix(x)
			if (!is.numeric(z)) { 
				stop("if x is a matrix, it must be a numeric matrix")
			}
			if (!missing(subset)) z <- z[subset, , drop = FALSE]
			if(robust) {
				if(!require(robustbase)) {
					stopifnot(require(MASS))
					covmat <- cov.mcd(x = z, cor = TRUE, 
							nsamp = "best", seed = seeds)
				}
				else covmat <- covMcd(x = z, cor = TRUE, 
							nsamp = "best", seed = seeds)
			}
			else covmat <- cov.wt(z, cor = TRUE)
        	}
		covmat$zz <- scale(z, TRUE, TRUE)
		covmat$n.obs <- nrow(z)
	}
# 	else stop("'covmat' is of unknown type")
	if(have.x) covmat$z <- z
	cormat <- cov2cor(covmat$cov)
	ev <- eigen(cormat, TRUE, TRUE)$values
	if(ev[ncol(cormat)] < sqrt(.Machine$double.eps)) {
		warning("sample covariance matrix is not positive definite\n",
			"continuing with an adjustment to make it so", call. = FALSE)
		covmat$cov <- FAiR_nearPD(covmat$cov, maxit = .Machine$integer.max)
		cormat <- cov2cor(covmat$cov)
	}
	covmat$cor <- cormat
	class(covmat) <- "list"
	return(covmat)
}

FAiR_PACE <- # communalities via brute force partitioning method; see Kano(1990)
function(Sigma, factors) { ## this takes forever in nontrivial problems
	temp.fun <- function(x) {
		the.list <- rep(list(x), factors)
		eg <- expand.grid(the.list)
		eg.sorted <- t(apply(eg, 1, FUN = function(z) {
			if(any(duplicated(z))) return(rep(NA_real_, length(x)))
			y <- x[!(x %in% z)]
			return(c(sort(z), sort(y)))
		}))
		eg <- unique(eg.sorted)
		eg <- eg[complete.cases(eg),1:factors,drop = FALSE]
		return(eg[1:(nrow(eg)/2),])
	}

	combos <- combn(1:nrow(Sigma), 2 * factors, FUN = temp.fun)

	dets <- apply(combos, 3, FUN = function(x) {
		uniques <- unique(c(x))
		apply(x, 1, FUN = function(two) {
			three <- uniques[!(uniques %in% two)]
			det(Sigma[two,three,drop = FALSE])
		})
	})
	communalities <- rep(NA_real_, nrow(Sigma))
	for(i in 1:nrow(Sigma)) {
		included <- apply(combos, 3, FUN = function(x) !(i %in% x))
		included_combos <- combos[,,included,drop = FALSE]
		included_dets <- dets[,included, drop = FALSE]
		best <- which.max(abs(included_dets))
		the.shelf <- ceiling(best / nrow(included_combos))
		the.row <- best %% nrow(included_combos)
		if(the.row == 0) the.row <- nrow(included_combos)
		uniques <- unique(c(included_combos[,,the.shelf]))
		two <- included_combos[the.row,,the.shelf]
		three <- uniques[!(uniques %in% two)]
		
		rho_21 <- Sigma[two, i]
		rho_31 <- Sigma[three,i]
		rho_23 <- Sigma[two,three]

		communalities[i] <- rho_31 %*% solve(rho_23) %*% rho_21
	}
	communalities <- ifelse(communalities >= 1, 1 - .Machine$double.eps, 
			 ifelse(communalities <= 0, 0 + .Machine$double.eps,
				communalities))
	names(communalities) <- colnames(Sigma)
	return(communalities)
}

FAiR_PACE_by_RGENOUD <- # this is much faster but not guaranteed to give the right answer
function(Sigma, factors, pp = paste(tempfile(), "PACE.txt", sep = ""), seeds, ...) {
	stopifnot(require(rgenoud))
	if(missing(seeds)) seeds <- c(formals(genoud)[c("unif.seed", "int.seed")])
	else if(length(seeds) == 1) seeds <- rep(seeds, 2)
	opt.fun <- function(par, first) {
		if(first %in% par)       return(c(0,0,0))
		if(any(duplicated(par))) return(c(1,0,0))
		return(c(1,1,abs(det(Sigma[     par[1:factors], 
						par[-c(1:factors)], drop = FALSE ]))))
	}
	communalities <- (0.5 * factors / ncol(Sigma)) / diag(solve(Sigma))
	for(i in 1:nrow(Sigma)) {
		proj.path <- paste(pp, i, sep = "_")
		if(i != 1) seeds <- c(formals(genoud)[c("unif.seed", "int.seed")])
		opt <- genoud(fn = opt.fun, nvars = 2 * factors, max = TRUE,
				Domains = cbind(rep(1.0, 2 * factors), nrow(Sigma)),
				boundary.enforcement = 2, lexical = TRUE, print.level = 0,
				data.type.int = TRUE, project.path = proj.path, 
				first = i, MemoryMatrix = FALSE, P3 = 0, 
				unif.seed = seeds[1], int.seed = seeds[2], ...)
		if(any(opt$value < sqrt(.Machine$double.eps))) {
			warning("submatrix almost singular, using SMC")
			next
		}

		two   <- opt$par[   1:factors,  drop = FALSE]
		three <- opt$par[-c(1:factors), drop = FALSE]

		rho_21 <- Sigma[two,  i]
		rho_31 <- Sigma[three,i]
		rho_23 <- Sigma[two,three, drop = FALSE]

		communalities[i] <- rho_31 %*% solve(rho_23) %*% rho_21	
	}
	communalities <- ifelse(communalities >= 1, 1 - sqrt(.Machine$double.eps), 
			 ifelse(communalities <= 0, 0 + sqrt(.Machine$double.eps),
				communalities))
	names(communalities) <- colnames(Sigma)
	return(communalities)
}

FAiR_make_starts_for_Rotate <- # makes pop.size starting values
function(A, pop.size) {
	n <- nrow(A)
	r <- ncol(A)
	A_swept <- sweep(A, 1, sqrt(rowSums(A^2)), FUN = "/", check.margin = FALSE)
	cutoff <- sort(combn(1:n, r, FUN = function(x) 
			det(crossprod(A_swept[x,]))), TRUE)[floor(pop.size / 2)]
	
	good <- combn(1:n, r, simplify = FALSE, 
		FUN = function(x) if(det(crossprod(A_swept[x,])) >= cutoff) 
					return(x) else return(NULL) )
	
	out <- array(NA, dim = c(r, r, 2 * floor(pop.size / 2)))
	count <- 0
	for(i in 1:length(good)) {
		if(is.null(good[[i]])) next
		count <- count + 1
		Tmat <- t(A_swept[good[[i]],])
		Tmat_svd <- svd(Tmat)
		out[,,count] <- Tmat
		out[,,pop.size - count + 1] <- Tmat_svd$u %*% t(Tmat_svd$v)
	}
	signs <- out[r,,]
	out <- matrix(c(out[-r,,]), ncol = pop.size)
	out <- rbind(out, signs)
	return(out)
}

FAiR_Landahl <- # (possibly) used to generate an orthogonal starting Tmat
function(r) {   # r is the number of factors
	k <- (r-1):1
	submatrix <- diag( sqrt(k / (k+1) ) )
	Tmat <- rbind( 1/sqrt(r), cbind( submatrix, -1 / sqrt( (k+1) * k) ) )
	for(k in 2:(r-1)) for(j in k:(r-1)) Tmat[k,j] <- Tmat[k,r]
	return(Tmat)
}

FAiR_make_Tmat <- # intelligently coerces a vector into an oblique transformation matrix
function(par) {
	factors <- sqrt(length(par))
	end <- factors^2 - factors

	# fill all but the last row of Tmat
	Tmat <- matrix(par[1:end], nrow = factors - 1, ncol = factors)
	colsums <- colSums(Tmat^2)
	Tmat <- rbind(Tmat, NA_real_)   # append another row to Tmat
	if(any(colsums > 1)) {          # decomplexify Tmat
		attributes(Tmat)$too.big <- max(colsums)
		colsums <- complex(length.out = factors, real = colsums)
		Tmat[factors,] <- sqrt(1 - colsums) * # get the signs for the last row
					ifelse(par[(end +1):length(par)] > 0, 1, -1)
		Tmat <- Re(Tmat)	# normalize the real part of the complex matrix
		Tmat <- sweep(Tmat, 2, sqrt(colSums(Tmat^2)), "/", FALSE)
	}
	else {
		attributes(Tmat)$too.big <- -1.0
		Tmat[factors,] <- sqrt(1 - colsums) * # get the signs for the last row
					ifelse(par[(end + 1):length(par)] > 0, 1, -1)
	}
	return(Tmat)
}

FAiR_Rotate <- 
function(par, A, criteria) {
	out <- rep(0, 1 + length(criteria))
	Tmat <- FAiR_make_Tmat(par)
	out[1] <- attributes(Tmat)$too.big
	if(out[1] != -1) return(out)
	Phi  <- crossprod(Tmat)
	out[2] <- criteria[[1]](Phi)
	if(out[2] != -1) return(out)

	rotmat_primary   <- Tmat %*% chol2inv(chol(Phi))
	rotmat_reference <- sweep(rotmat_primary, 2, sqrt(colSums(rotmat_primary^2)), 
				FUN = "/", check.margin = FALSE)

	Upsilon <- A %*% rotmat_reference
	
	# calculate all criteria
	for(i in 2:length(criteria)) {
		environment(criteria[[i]]) <- environment()
		out[i + 1] <- criteria[[i]]()
		if(out[i + 1] != -1) return(out)
	}
}

# FAiR_extended_vector_plot <- # rethink this
# function(x, which.cols = 2:4, ...) {
# 	## handle cases where length(which.cols) > 4 or ncol(x) <= 2
# 	x <- sweep(x, 1, x[,1], FUN = "/")
# 	if(ncol(x) == 3 | length(which.cols) == 2) {
# 		par(pty = "s")
# 		plot(x[,which.cols[1:2]])
# 	}
# 	else if(ncol(x) == 4 | length(which.cols) == 3) {
# 	        stopifnot(require(rgl))
# 		the.range <- range(x[,which.cols[1:3]])
# 		plot3d(x[,which.cols[1:3]], xlim = the.range, ylim = the.range, 
# 			zlim = the.range, type = "s",
# 			radius = .05) # Do xlab, ylab, and zlab
# 	}
# 	else stop("this should not happen, please report bug to maintainer")
# }

FAiR_Tmat_approx <-
function(Tmat) {
	Tmat_svd <- svd(Tmat)
	return(Tmat_svd$u %*% t(Tmat_svd$v))
}


FAiR_fix_coefficients <- # this function is for setting exact zeros in a SEFA
function(coefs, cormat = diag(ncol(coefs)), zeros = ncol(coefs), quasi_Yates = FALSE,
	weak_Thurstone = FALSE, Butler = FALSE, row_complexity = NA) {

	D_inv <- sqrt(rowSums(backsolve(chol(cormat), diag(ncol(cormat)))^2))
	Upsilon <- sweep(coefs, 2, D_inv, FUN = "/", check.margin = FALSE)

	if(!is.na(row_complexity[1])) {
		if(length(row_complexity) == 1) {
			Upsilon <- t(apply(abs(Upsilon), 1, FUN = function(x) {
				ranks <- order(x, decreasing = TRUE)
				x[ranks > row_complexity] <- 0
				return(x)
			}))
		}
		else {
			joined <- cbind(row_complexity, abs(Upsilon))
			Upsilon <- t(apply(joined, 1, FUN = function(x) {
				row_complexity <- x[1]
				x <- x[-1]
				ranks <- order(x, decreasing = TRUE)
				x[ranks > row_complexity] <- 0
				return(x)
			}))
		}
		coefs <- sweep(Upsilon, 2, D_inv, FUN = "*", check.margin = FALSE)
	}

	else if(weak_Thurstone) { # r zeros per factor, all of complexity r - 1
		Upsilon.abs <- abs(Upsilon)
		orderer <- order(c(Upsilon.abs))
		changed <- rep(0, ncol(coefs))
		for(i in 1:length(orderer)) {
			dims <- which(coefs == coefs[orderer[i]], arr.ind = TRUE)
			if(changed[dims[1,2]] == ncol(coefs)) next
			if(any(!coefs[dims[1,1], ])) next
			coefs[orderer[i]] <- 0
			changed[dims[1,2]] <- changed[dims[1,2]] + 1
			if(all(changed == ncol(coefs))) break
		}

		if(all(zeros <= ncol(coefs))) return(coefs) # can add more zeros
		Upsilon <- sweep(coefs, 2, D_inv, "/", check.margin = FALSE)
	}

	else if(Butler) { ## RETHINK THIS AGAIN
		# First make distinguishability weights following Yates (1987)
		FC <- coefs * (coefs %*% cormat)
		w <- colSums(FC) / sum(FC) # similar to equation 118 in Yates
		normalizer <- sqrt(rowSums(FC))
		cosines <- apply(coefs, 1, FUN = "%*%", w) / normalizer # 119 in Yates
		factors <- ncol(coefs) # dw is equation 78a in Yates
		dw <- pmin(1, cosines^2 * (1 - cosines^2) * factors^2 / (factors - 1) )

		# Then unifactorialize the test with the biggest dw for each factor
		biggest <- apply(FC, 1, which.max)
		ordered <- order(dw, decreasing = TRUE)
		unifactorial <- rep(FALSE, factors)
		for(i in 1:nrow(coefs)) {
			mark <- biggest[ordered[i]]
			if(!unifactorial[mark]) {
				coefs[ordered[i], -mark] <- 0
				coefs[ordered[i],  mark] <- abs(coefs[ordered[i], mark])
				unifactorial[mark] <- TRUE
			}
			if(all(unifactorial)) break
		}

		if(all(zeros < ncol(coefs))) return(coefs) # can add more zeros
		Upsilon <- sweep(coefs, 2, D_inv, "/", check.margin = FALSE)

# 		MOST RECENT WAY STARTS NOW
# 		FC <- coefs * (coefs %*% cormat)
# 		maxs <- apply(FC, 2, max)
# 		for(i in 1:ncol(coefs)) {
# 			coefs[maxs[i], -i] <- 0
# 			coefs[maxs[i],  i] <- abs(coefs[maxs[i], i])
# 		}

# 		OLD WAY STARTS NOW
# 		FC <- coefs * (coefs %*% cormat)
# 		ratios <- sweep(FC, 1, rowSums(FC), FUN = "/")
# 		maxs <- apply(ratios, 2, FUN = which.max)
# 		for(i in 1:ncol(coefs)) { 
# 			coefs[maxs[i], -i] <- 0
# 			coefs[maxs[i],  i] <- abs(coefs[maxs[i],  i])
# 		}

# 		OLD, OLD WAY STARTS NOW
# 		D_inv <- sqrt(diag(chol2inv(chol(cormat))))
# 		Upsilon.abs <- abs(sweep(coefs, 2, D_inv, FUN = "/", FALSE))
# 		rowsums <- rowSums(Upsilon.abs)
# 		distinguished <- apply(sweep(-Upsilon.abs, 2, rowsums, FUN = "+", FALSE),
# 					2, FUN = which.min)
# 		for(i in 1:ncol(coefs)) coefs[distinguished[i], -i] <- 0
# 		if(all(zeros < ncol(coefs))) return(coefs)
	}

	else if(quasi_Yates) { # cohyperplanarity without minimal collinearity
		FC <- as.matrix( coefs * (coefs %*% cormat) )
		for(i in 1:ncol(coefs)) {
			diffs <- sweep(FC, 1, FC[,i], FUN = "-", check.margin = FALSE)
			coefs[apply(diffs, 2, which.max)[-i], i] <- 0
		}
		if(all(zeros < ncol(coefs))) return(coefs) # can add more zeros
		Upsilon <- sweep(coefs, 2, D_inv, "/", check.margin = FALSE)
	}

	if(length(zeros) == 1) zeros <- rep(zeros, ncol(coefs))
	Upsilon <- abs(Upsilon)
	for(i in 1:ncol(coefs)) {
		ordered <- order(Upsilon[,i])
		coefs[ordered[1:zeros[i]],i] <- 0
	}
	return(coefs) # a matrix
}

FAiR_check_coefficients <- # this function does the submatrix rank check for Howe (1955)
function(coefs, threshold) {
	not.duds <- apply(coefs, 1, FUN = any)
	for(i in 1:ncol(coefs)) {
		subcoefs <- coefs[coefs[,i] == 0 & not.duds, ,drop = FALSE]
		if(nrow(subcoefs) == 0)         return(-.Machine$double.xmax)
		if(any(is.infinite(subcoefs)))  return(-.Machine$double.xmax)
		if( sum(svd(subcoefs)$d > threshold) != (ncol(subcoefs) - 1) ) {
						return(-.Machine$double.xmax)
		}
	}
	return(1.0)
}

FAiR_opt2FAobject <- 
function(opt, FAobject, seeds) {
	manifest <- FAobject@manifest
	S <- manifest$cor
	Lambda <- coef(FAobject)

	Tmat <- FAiR_make_Tmat(opt$par)

	Phi  <- crossprod(Tmat)
	rotmat_primary   <- Tmat %*% chol2inv(chol(Phi))
	rotmat_reference <- sweep(rotmat_primary, 2, sqrt(colSums(rotmat_primary^2)), 
				FUN = "/", check.margin = FALSE)

	Phi <- crossprod(Tmat)
	Psi <- crossprod(rotmat_reference)
	D   <- diag(diag(crossprod(Tmat, rotmat_reference)))

	PP <- Lambda %*% rotmat_primary
	RS <- Upsilon <- Lambda %*% rotmat_reference
	PS <- Lambda %*% Tmat
	RP <- Lambda %*% rotmat_reference %*% chol2inv(chol(Psi))
	FC <- PP * PS
	sorter <- order(colSums(FC), decreasing = TRUE)

	loadings <- array(cbind(PP, RS, PS, RP, FC), dim = c(dim(PP), 5),
				dimnames = list(rownames(Lambda), NULL, 
				c("PP", "RS", "PS", "RP", "FC")))

        correlations <- array(cbind(Phi, Psi, D), dim = c(rep(ncol(D), 2), 3),
                                dimnames = list(NULL, NULL, c("PF", "RF", "PR")))

	signs <- sign(colSums(loadings[,,1]))
	signs[signs == 0] <- 1
	if(any(signs != 1)) {
		loadings <- sweep(loadings, 2, signs, FUN = "*")
		loadings[,,5] <- FC
	        correlations  <- sweep(sweep(correlations, 1, signs, FUN = "*"),
         	                                           2, signs, FUN = "*")
		Phi <- sweep(sweep(Phi, 1, signs, FUN = "*"), 2, signs, FUN = "*")
	}

	par <- sweep(loadings[,,1], 1, sqrt(diag(manifest$cov)), FUN = "*", FALSE)
	par <- c(Phi[upper.tri(Phi)], par, FAobject@uniquenesses * diag(manifest$cov))

	if(is.numeric(manifest$n.obs) && FAobject@method == "MLE") {
		num_free <- length(loadings[,,1])
		mark_start <- 1
		mark_end <- 0.5 * ncol(Lambda) * (ncol(Lambda) - 1)
		p <- nrow(S)
		beta_list <- list(beta = loadings[,,1], free = rep(TRUE, num_free),
					num_free = num_free,
					select = c(rep(FALSE, mark_end), rep(TRUE,  num_free),
									 rep(FALSE, p)) )
		Theta2_list <- list(Theta2 = diag(nrow(S)), 
					select = c(rep(FALSE, mark_end + num_free), rep(TRUE, p)))
		new_restrictions <- new("restrictions.1storder", factors = c(ncol(Lambda),0),
					nvars = length(par), dof = FAobject@restrictions@dof,
					Domains = cbind(par, par), model = "CFA",
					method = "MLE", Phi = diag(rep(0.5, ncol(Lambda))),
					beta = beta_list, Theta2 = Theta2_list,
					criteria = list(FAiR_criterion_llik))
		vcov   <- FAiR_BHHH(par, new_restrictions, manifest)		
		zstats  <- par / sqrt(diag(vcov))
		zstats_Phi <- 0 * Phi
		zstats_Phi[upper.tri(Phi)] <- zstats[mark_start:mark_end]
		zstats_Phi <- zstats_Phi + t(zstats_Phi)
		diag(zstats_Phi) <- NA_real_
		zstats_Phi <- zstats_Phi[sorter,sorter]
		zstats_beta <- Lambda * NA_real_
		mark_start <- mark_end + 1
		mark_end <- mark_end + num_free
		fill <- new_restrictions@beta$free
		zstats_beta[fill] <- zstats[mark_start:mark_end]
		zstats_beta <- zstats_beta[,sorter]
		mark_start <- mark_end + 1
		mark_end <- length(zstats)
		zstats_Theta2 <- zstats[mark_start:mark_end]
		zstats <- list(Phi = zstats_Phi, beta = zstats_beta, 
				Theta2 = zstats_Theta2)
	}
	else {
		warning("z-statistics could not be calculated because the number of observations",
			" is unknown")
		vcov <- matrix(NA_real_, nrow = length(par), ncol = length(par))
		zstats <- list()
	}

	FAobject@loadings <- loadings[,sorter,]
	FAobject@correlations <- correlations[sorter,sorter,]
	attributes(FAobject@correlations)$orthogonal <- FALSE

	trans_mats <- array(cbind(rotmat_primary, rotmat_reference), 
				dim = c(rep(ncol(D), 2), 2),
				dimnames = list(NULL, NULL, c("primary", "reference")))

	FAobject@trans_mats <- trans_mats[sorter,sorter,]
	colnames(FAobject@loadings) <-  colnames(FAobject@correlations) <-
		rownames(FAobject@correlations) <- colnames(FAobject@trans_mats) <- 
		colnames(FAobject@trans_mats) <- paste("F", 1:ncol(D), sep = "")

	scores <- FAobject@call$scores
	if(is.null(scores)) scores <- "none"
 	FAobject@scores  <- FAiR_scores(scores, manifest$zz, beta, Phi, 
					FAobject@uniquenesses, S)
	FAobject@rgenoud$transformation <- opt
	FAobject@seeds   <- rbind(extraction     = c(FAobject@seeds),
				  transformation = seeds)
	return(FAobject)
}

FAiR_BHHH <-
function(par, restrictions, manifest) {
	free <- par != 0
	free_sum <- sum(free)
	if(restrictions@method != "MLE") {
		warning("BHHH can only be calculated under maximum likelihood")
		return(matrix(NA_real_, free_sum, free_sum))
	}
	if(manifest$n.obs < sum(free)) {
		warning("BHHH can only be calculated when the number of observations",
			"is greater than or equal to the number of free parameters")
		return(matrix(NA_real_, free_sum, free_sum))
	}
	if(is.null(manifest$zz)) {
		warning("BHHH cannot be calculated because the raw data were not",
			" passed to Factanal()")
		return(matrix(NA_real_, free_sum, free_sum))
	}
	restrictions@model <- "CFA"
	z <- sweep(manifest$zz, 2, sqrt(diag(manifest$cov)), FUN = "*")
	helper2 <- bfgs_helpS4(par, object = restrictions, done = TRUE, 
			       S = manifest$cor, lower = sqrt(.Machine$double.eps))
	G <- t(apply(z, 1, FUN = function(x) {
			gr_fitS4(par, object = restrictions, helper2, 
				  S = tcrossprod(x), lower = sqrt(.Machine$double.eps))
		}))
	G <- G[,free] * sqrt(0.5)
	GG_inv <- try(chol2inv(chol(crossprod(G))))
	if(is.matrix(GG_inv)) return(GG_inv)
	else {
		warning("variance-covariance appears not to be positive definite")
		GG_inv <- try(solve(crossprod(G)))
		if(is.matrix(GG_inv)) return(GG_inv)
		else return(matrix(NA_real_, free_sum, free_sum))
	}
}

FAiR_gradient <- # calculates the gradient of the log-likelihood wrt beta, Theta2, and Phi
function(cormat, covmat = NULL, do_Phi = TRUE, Phi, beta, Theta2) {
	if(!is.null(covmat)) {
		S <- covmat$cov
		scale <- (covmat$n.obs - 1) / 2
	}
	else {
		S <- cormat
		scale <- 1
	}

	C <- crossprod(chol(Phi) %*% t(beta)) + Theta2
	C_inverse <- chol2inv(chol(C))

	middle <- C_inverse %*% (S - C) %*% C_inverse

	dF_dbeta <- middle %*% beta %*% Phi * scale * 2
	dF_dTheta2 <- diag(middle) * scale

	out <- list(dF_dbeta = dF_dbeta, dF_dTheta2 = dF_dTheta2)
	if(do_Phi) out$dF_dPhi <- t(beta) %*% middle %*% beta * scale
	return(out)
}

FAiR_beta_under_orthogonality <-
function(x, fixed, communalities) {
	beta <- fixed
	beta[is.na(fixed)] <- x
	beta <- sweep(beta, 1, sqrt(rowSums(beta^2) / communalities), 
			FUN = "/", check.margin = FALSE)
	signs <- sign(colSums(beta))
	beta <- sweep(beta, 2, signs, FUN = "*", check.margin = FALSE)
	return(beta[is.na(fixed)])
}

FAiR_lexical_driver <-
function(par, object, S, Delta = NULL, Xi = NULL, lower) {
	out <- rep(0, 1 + length(object@criteria))
	Phi <- object@Phi
	beta <- object@beta$beta
	Theta2 <- object@Theta2$Theta2
	C <- crossprod(chol(Phi) %*% t(beta)) + Theta2

	ev <- eigen(C, TRUE, TRUE)$values
	if(ev[ncol(C)] < lower) {
		out[1] <- ev[ncol(C)]
		return(out)
	}
	else out[1] <- 1.0

	marker <- 1L

	# calculate all remaining criteria
	for(i in 1:length(object@criteria)) {
		environment(object@criteria[[i]]) <- environment()
		if( (out[marker + i] <- object@criteria[[i]]()) != 1 ) return(out)
	}
}

FAiR_nearPD <- # makes x into a positive definite correlation matrix
function(x, corr = TRUE, eig.tol = 1e-06, conv.tol = 1e-07, 
	posd.tol = sqrt(.Machine$double.eps), do2eigen = TRUE, maxit = 100) {
	## This function is slightly modified from Matrix:::nearPD, which is
	## Copyright (2007) Jens Oehlschl�gel and is licensed under the GPL V2+
	e <- eigen(x, symmetric = TRUE, only.values = TRUE)$values
	n <- ncol(x)
	if(e[n] > posd.tol * e[1]) {
		attributes(x)$ev <- 1.0
		return(x)
	}
	U <- x
	U[] <- 0
	X <- x
	iter <- 0
	converged <- FALSE
	conv <- Inf
	inorm <- function(x) max(rowSums(abs(x)))
	while (iter < maxit) {
		Y <- X
		T <- Y - U
		e <- eigen(Y, symmetric = TRUE)
		Q <- e$vectors
		d <- e$values
		D <- diag(d)
		p <- d > eig.tol * d[1]
		X <-      Q[, p, drop = FALSE] %*% D[p, p, drop = FALSE] %*% 
			t(Q[, p, drop = FALSE])
		U <- X - T
		X <- (X + t(X))/2
		if (corr) diag(X) <- 1
		conv <- inorm(Y - X) / inorm(Y)
		converged <- (conv <= conv.tol)
		if(converged) break
		iter <- iter + 1
	}
	X <- (X + t(X))/2
	if (do2eigen) {
		e <- eigen(X, symmetric = TRUE)
		d <- e$values
		Eps <- posd.tol * abs(d[1])
		if (d[n] < Eps) {
			d[d < Eps] <- Eps
			Q <- e$vectors
			o.diag <- diag(X)
			X <- Q %*% (d * t(Q))
			D <- sqrt(pmax(Eps, o.diag)/diag(X))
			X[] <- D * X * rep(D, each = n)
		}
	}
	if (corr) X <- cov2cor(X)
	attributes(X)$ev <- min(d)
	return(X)
}

FAiR_stats <- 
function(object, conf.level = 0.9) {
	MLE <- object@restrictions@method == "MLE" 
	N  <- n.obs <- object@manifest$n.obs
	S  <- model.matrix(object)
	n  <- nrow(S)
	C  <- fitted(object) + diag(object@uniquenesses)
	if(MLE) {
		llik  <- logLik(object)
		df <- attributes(llik)$df
		llik_saturated <- -0.5 * (n.obs - 1) * (log(det(S)) + nrow(S))
		df_0 <- 0.5 * n * (n + 1)
		p.value <- pchisq(2 * deviance(object), df_0 - df, lower.tail = FALSE)
		BIC_saturated <- -2 * llik_saturated + df_0 * log(n.obs)
		sds <- sqrt(diag(object@manifest$cov))
		C_unscaled <- sweep(sweep(C, 1, sds, FUN = "*"),
					     2, sds, FUN = "*")
		llik_unscaled <- -log(det(C_unscaled)) - 
				  sum(diag(object@manifest$cov %*% solve(C_unscaled)))
		if(all(!is.na(vcov(object)))) { 
			info <- try(chol2inv(chol(vcov(object))))
			if(!is.matrix(info)) info <- try(solve(vcov(object)))
			if(!is.matrix(info)) info <- NA_real_ * vcov(object)
		}
		else info <- NA_real_ * vcov(object)
		SIC <- -llik_unscaled * (N - 1) + n * log(N) + log(det(info))
		SIC <- SIC / 2
	}
	else {
		llik <- llik_saturated <- BIC_saturated <- SIC <- p.value <- NA_real_
		RMSEA <- c(rep(NA_real_, 3), conf.level)
	}

	## The following is slightly modified from sem:::summary.sem, which is
	## Copyright 2007 John Fox and is licensed under the GPL V2+
	invC <- chol2inv(chol(C))
	CNull <- diag(diag(S))
	invCNull <- diag(1/diag(S))
	chisqNull <- log(det(CNull)) - log(det(S)) + 
			crossprod(as.vector(S), as.vector(invCNull)) - n
	chisqNull <- chisqNull * (N - 1)
	chisq <- log(det(C)) - log(det(S)) + 
			crossprod(as.vector(S), as.vector(invC)) - n
	chisq <- chisq * (N - 1)
	dfNull <- n*(n - 1)/2
	df <- df.residual(object)
	if(!MLE) chisq <- NA_real_
	CSC <- invC %*% (S - C)
	CSC <- CSC %*% CSC
	CS  <- invC %*% S
	CS  <- CS %*% CS
	GFI  <- 1 - sum(diag(CSC)) / sum(diag(CS))
        AGFI <- 1 - (n * (n + 1) / (2 * df)) * (1 - GFI)
        NFI  <- (chisqNull - chisq) / chisqNull
        NNFI <- (chisqNull / dfNull - chisq / df) / (chisqNull / dfNull - 1)
        L1 <- max(chisq - df, 0)
        L0 <- max(L1, chisqNull - dfNull)
        CFI <- 1 - L1/L0
        RMSEA <- sqrt(max(chisq/(df * (N - 1)) - 1/(N - 1), 0))
        tail <- (1 - conf.level)/2 
        max <- n.obs

        if(MLE) {
		while (max > 1){
			res <- optimize(function(lam) (tail - pchisq(chisq, df,
				ncp=lam))^2, interval=c(0, max))
			if (sqrt(res$objective) < tail/100) break
			max <- max/2
		}
        	lam.U <- if (max <= 1) NA_real_ else res$minimum
		max <- max(max, 1)

        	while (max > 1){
			res <- optimize(function(lam) (1 - tail - pchisq(chisq, df,
					ncp=lam))^2, interval=c(0, max))
			if (sqrt(res$objective) < tail/100) break
			max <- max/2
		}
		lam.L <- if (max <= 1) NA_real_ else res$minimum
		RMSEA.U <- sqrt(lam.U / ((N - 1)*df))
		RMSEA.L <- sqrt(lam.L / ((N - 1)*df))
		RMSEA <- c(RMSEA, RMSEA.L, RMSEA.U, conf.level)
	}

	return( list(llik = llik, p.value = p.value,
			GFI = GFI, AGFI = AGFI, NFI = NFI, NNFI = NNFI, CFI = CFI,
			RMSEA = RMSEA, BIC = BIC(object), 
			BIC_saturated = BIC_saturated, SIC = SIC) )
}

FAiR_get_loadings <-
function(Psi, S, q) {
	sc <- diag(1/sqrt(Psi))
	Sstar <- sc %*% S %*% sc
	E <- eigen(Sstar, symmetric = TRUE)
	L <- E$vectors[, 1:q, drop = FALSE]
	load <- L %*% diag(sqrt(pmax(E$values[1:q] - 1, 0)), q)
	diag(sqrt(Psi)) %*% load
}

FAiR_triads <- # estimates the coefficents for a 3x3 correlation matrix 
function(Phi) {
	cors <- Phi[lower.tri(Phi)]
	triads <- prod(cors) / cors^2
	Lambda <- (sqrt(triads) * sign(cors))[c(2,3,1)]
	return(Lambda)
}

FAiR_tetrads <- # estimates the coefficents for a 4x4 correlation matrix 
function(Phi) {
	A <- matrix(NA, nrow = 4, ncol = 3)
	A[1,] <- c(Phi[1,2] * Phi[1,3] / Phi[2,3],
		   Phi[1,2] * Phi[1,4] / Phi[2,4],
		   Phi[1,3] * Phi[1,4] / Phi[3,4])
	A[2,] <- c(Phi[2,1] * Phi[2,3] / Phi[1,3],
		   Phi[2,1] * Phi[2,4] / Phi[1,4],
		   Phi[2,3] * Phi[2,4] / Phi[3,4])
	A[3,] <- c(Phi[3,1] * Phi[3,2] / Phi[1,2],
		   Phi[3,1] * Phi[3,4] / Phi[1,4],
		   Phi[3,2] * Phi[3,4] / Phi[2,4])
	A[4,] <- c(Phi[4,1] * Phi[4,2] / Phi[1,2],
		   Phi[4,1] * Phi[4,3] / Phi[1,3],
		   Phi[4,2] * Phi[4,3] / Phi[2,3])
	## What do I do with A? See which column fits best?
	return(A)
}

FAiR_numeric_gradient <- 
function(par, object, helper, S) {
	gradient.env <- new.env()
	assign("par", par, envir = gradient.env)
	assign("object", object, envir = gradient.env)
	assign("helper", helper, envir = gradient.env)
	assign("S", S, envir = gradient.env)
	gradient <- as.real(attr(numericDeriv(quote(bfgs_fitS4(par, object,
			helper, S, lower = sqrt(.Machine$double.eps))), 
			theta = c("par"), gradient.env), "gradient"))
	return(gradient)
}

FAiR_informatrix <-
function(par, object, helper, S) {
	A <- dC_dpar %*% kronecker(C_inv, C_inv) %*% dC_dpar
	
}

FAiR_Browne1974 <-
function(C) {
	vec  <- as.vector(C)
	vecs <- C[lower.tri(C, diag = TRUE)]

	Kp <- matrix(0, nrow = length(vec), ncol = length(vecs))
	mark <- 1
	for(i in 1:length(vec)) {
		colmark <- which(vecs == vec[i])
		if(vec[i] == 1) {
			colmark <- colmark[mark]
			mark <- mark + 1
		}
		Kp[i,colmark] <- 1
	}
	Kp_prime <- solve(crossprod(Kp)) %*% t(Kp)
	stopifnot(all.equal(vec, c(Kp %*% vecs), check.attributes = FALSE))
	
	Gamma_N <- 2 * Kp_prime %*% kronecker(C, C) %*% t(Kp_prime)
}

FAiR_EM <-
function() {
# 	stopifnot(require(Design))
	stop("write this function")

	notconverged <- TRUE
	while(notconverged) { # break if a constraint is not met
		## E step
	
		## Get unconstrained beta
	
		## Transform to constrained beta
	
		## Get Phi or Delta and Xi
	
		## Get uniquenesses

		## Check convergence

		## Move temp stuff to real stuff
	}
}

FAiR_Hessian <-
function(par, object, S, lower) {
	## This function is slightly modified from stats:::optim, which is 
        ## Copyright 1995-2007 R Core Development Team and licensed under GPL V2+

	object@model <- "CFA"
	help <- bfgs_helpS4(par, object, done = TRUE, S, lower)
	fn1 <- function(par, helper = help) bfgs_fitS4(par, object, helper, S, lower)
	gr1 <- function(par, helper = help)   gr_fitS4(par, object, helper, S, lower)
	if(object@method != "MLE") gr1 <- NULL

	con <- list(trace = 0, fnscale = -1, parscale = rep.int(1, 
		length(par)), ndeps = rep.int(0.001, length(par)), maxit = 100, 
		abstol = -Inf, reltol = sqrt(.Machine$double.eps), alpha = 1, 
		beta = 0.5, gamma = 2, REPORT = 10, type = 1, lmm = 5, 
		factr = 1e+07, pgtol = 0, tmax = 10, temp = 10)

	hess <- .Internal(optimhess(par, fn1, gr1, con))
	hess <- 0.5 * (hess + t(hess))
	return(hess)
}

FAiR_scores <-
function(scores, zz, beta, Phi, uniquenesses, S) {
	if(scores == "none") return(matrix(NA_real_, 1, 1))
	else if(is.null(zz)) {
		warning("could not calculate factor scores because the raw data on the",
			" outcome variables was not supplied")
		return(matrix(NA_real_, 1, 1))
	}
	Theta2inv <- diag(1/uniquenesses)
	S_inv <- chol2inv(chol(S))
	if(scores == "regression") {
		B <- solve(S, beta) %*% Phi
	}
	else if(scores == "Bartlett") { # check if same
		B <- Theta2inv %*% beta %*% solve(t(beta) %*% Theta2inv %*% beta)
	}
	else if(scores == "Thurstone") {
		B <- S_inv %*% beta %*% Phi
	}
	else if(scores == "Ledermann") {
		B <- Theta2inv %*% beta %*% solve( t(beta) %*% Theta2inv %*% beta +
			chol2inv(chol(Phi)) )
	}
	else if(scores == "Anderson-Rubin") {
		B <- Theta2inv %*% beta %*% diag(diag( t(beta) %*% Theta2inv %*% S %*% 
						Theta2inv %*% beta )^(-1/2))
	}
	else if(scores == "McDonald") {
		N <- t(chol(Phi)) # check transposition
		B <- Theta2inv %*% beta %*% N %*% diag(diag(t(N) %*% t(beta) %*% 
			Theta2inv %*% S %*% Theta2inv %*% beta %*% N)^(-1/2)) %*% t(N)
	}
	else if(scores == "Krinjen") {
		Phi_half <- diag(diag(Phi)^(1/2))
		B <- S_inv %*% beta %*% Phi_half %*% diag(diag(Phi_half %*%
			t(beta) %*% S_inv %*% beta %*% Phi_half)^(-1/2)) %*% Phi_half
	}
	else if(scores == "Takeuchi") {
		B <- S_inv %*% beta %*% diag(diag(t(beta) %*% S_inv %*% beta)^(-1/2))
	}
	else if(scores == "Harman") {
		B <- beta %*% chol2inv(chol(crossprod(beta)))
	}
	sc <- zz %*% B
	return(sc)
}

FAiR_indeterminator <-
function(model, factors, zeros1, zeros2, nonzeros1, nonzeros2) {
	if(model == "SEFA") {
		if(any(zeros1 < (factors[1] - 1) )) return("indeterminate1")
		else if(factors[2] == 0) {
			if(any(zeros1 >= factors[1]) | nonzeros1 > 0) return("determined")
			else return("minimal")
		}
		else if(any(zeros2 < (factors[2] - 1) )) return("indeterminate2")
		else return("determined")
	}
	else { # CFA
		if(any(zeros1 < (factors[1] - 1) )) return("indeterminate1")
		else if(factors[2] == 0) return("determined")
		else if(any(zeros2 < (factors[2] - 1) )) return("indeterminate2")
		else return("determined")
	}
}

FAiR_show_coef <- 
function(object, level) {
	hc <- 1.5
	if(level == 1) {
		fix_args <- object@beta$fix_beta_args
		coefs <- object@beta$beta
		select <- object@beta$select
	}
	else {
		fix_args <- object@Delta$fix_Delta_args
		coefs <- object@Delta$Delta
		select <- object@Delta$select
	}
	any_fixed <- apply(coefs, 1, FUN = function(x) any(!is.na(x)))
	if(any(any_fixed)) {
		cat("Fixed coefficients at level ", level, "\n")
		print(coefs[any_fixed,,drop = FALSE])
	}
	else cat("No fixed coefficients at level ", level, "\n")

	Domains_coefs <- object@Domains[select,,drop = FALSE]
	any_bounds <- apply(Domains_coefs, 1, FUN = 
			function(x) {
				out <- (x[1] != -hc) | (x[2] != hc)
				return(out)
			})

	if(any(any_bounds)) {
		cat("\nBounded coefficients at level ", level, "\n")
		print(Domains_coefs[any_bounds,,drop = FALSE])
	}
	else cat("\nAll coefficients at level ", level, "in [", -hc, ",", hc, "]\n")
	return(NULL)
}

FAiR_show_sefa <-
function(coef_args, level) {
	cat("\nZeros per factor at level ", level, "\n")
	mat <- matrix(coef_args$zeros, nrow = 1)
	rownames(mat) <- "zeros"
	colnames(mat) <- paste("F", 1:ncol(mat), sep = "")
	print(mat)

	if(is.na(coef_args$row_complexity[1])) {
		if(length(coef_args$row_complexity) == 1) {
			row_complexity <- coef_args$row_complexity
			if(is.na(row_complexity)) row_complexity <- ncol(mat)
			cat("All outcomes are of maximum complexity ", 
			    row_complexity, "\n")
		}
		else cat("Outcomes have various complexities\n")
	}
	else {
		cat("\nStipulations on zeros at level ", level, "\n", 
		if(coef_args$quasi_Yates) "\tEncourage cohyperplanarity\n" else NULL,
		if(coef_args$weak_Thurstone)"\tWeak simple structure\n" else NULL,
		if(coef_args$Butler) "\tUnifactorial basis\n" else NULL, "\n")
	}
	return(NULL)
}

FAiR_show_constraints <-
function(criteria) {
	items <-  c(
		paste("\tReference factors at level 2 must have more effective",
		      "variane than do the primary factors at level 1\n"),
		paste("\tPrimary factors at level 2 must have more effective",
		      "variane than do the primary factors at level 1\n"),
		paste("\tReference factors at level 1 must have more effective",
		      "variance than does the battery as a whole\n"),
		paste("\tPrimary factors at level 1 must have more effective",
		      "variane than does the battery as a whole\n"),
		"\tNo suppressor variables at level 2\n",
		"\tNo suppressor variables at level 1\n",
		paste("\tReference factors at level 2 must have more generalized",
		      "variance than do the primary factors at level 2\n"),
		paste("\tReference factors at level 1 must have more generalized",
		      "variance than do the primary factors at level 1\n"),
		paste("\tTests in hyperplanes have more effective variance than",
		      "does the battery as a whole\n"),
		"\tlog-likelihood\n", "\tYates (1987) weighted least squares\n")
	names(items) <- c("evRF_2nd", "evPF_2nd", "evRF_1st", "evPF_1st", 
			  "no_suppressors_2nd", "no_suppressors_1st", 
			  "dets_2nd", "dets_1st", "cohyperplanarity",
			  "llik", "YWLS")
	items <- items[names(items) %in% names(criteria)]
	names(items) <- NULL
	cat("\nLexical criteria:\n")
	cat(items)
	return(NULL)
}
