#     This file is part of FAiR, a program to conduct Factor Analysis in R
#     Copyright 2008 Benjamin King Goodrich
#     Some portions of this code are Copyright 1995-2007 R Core Development Team and
#     Copyright 2007 John Fox, are licensed under the GPL version 2 or later, and are
#     indicated so in the respective functions
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


## This file defines S4 methods

## NOTE: This file is meant to be read with 90 columns with 8 space tabs

## These six are S4 generics but probably would be of limited use outside FAiR
## To add a model to FAiR, each of the following are needed.
setGeneric("fitS4", def = function(par, object, ...)
		standardGeneric("fitS4"),       useAsDefault = FALSE, package = "FAiR"
# 	showMethods("fitS4")
# 	Function: fitS4 (package FAiR)
# 	par="numeric", object="restrictions.1storder"
# 	par="numeric", object="restrictions.2ndorder"
# 	par="numeric", object="restrictions.factanal"
# 	par="numeric", object="restrictions.general"
# 	par="numeric", object="restrictions.orthonormal"
)
setGeneric("gr_fitS4", def = function(par, object, helper, ...) 
		standardGeneric("gr_fitS4"),    useAsDefault = FALSE, package = "FAiR"
# 	showMethods("gr_fitS4")
# 	Function: gr_fitS4 (package FAiR)
# 	par="numeric", object="restrictions.1storder", helper="list"
# 	par="numeric", object="restrictions.2ndorder", helper="list"
# 	par="numeric", object="restrictions.factanal", helper="ANY"
# 	par="numeric", object="restrictions.general", helper="list"
# 	par="numeric", object="restrictions.orthonormal", helper="ANY"
)
setGeneric("bfgs_fitS4",  def = function(par, object, helper, ...)
                standardGeneric("bfgs_fitS4"),  useAsDefault = FALSE, package = "FAiR"
# 	showMethods("bfgs_fitS4")
# 	Function: bfgs_fitS4 (package FAiR)
# 	par="numeric", object="restrictions.1storder", helper="list"
# 	par="numeric", object="restrictions.2ndorder", helper="list"
# 	par="numeric", object="restrictions.factanal", helper="ANY"
# 	par="numeric", object="restrictions.general", helper="list"
# 	par="numeric", object="restrictions.orthonormal", helper="ANY"
)
setGeneric("bfgs_helpS4", def = function(initial, object, done = FALSE, ...)
                standardGeneric("bfgs_helpS4"), useAsDefault = FALSE, package = "FAiR",
# 	showMethods("bfgs_helpS4")
# 	Function: bfgs_helpS4 (package FAiR)
# 	initial="numeric", object="restrictions.1storder", done="logical"
# 	initial="numeric", object="restrictions.2ndorder", done="logical"
# 	initial="numeric", object="restrictions.general", done="logical"
# 	initial="numeric", object="restrictions", done="logical"
# 	initial="numeric", object="restrictions.orthonormal", done="logical"
)
setGeneric("create_FAobject", def = function(restrictions, ...) standardGeneric("create_FAobject"),
		useAsDefault = FALSE, package = "FAiR"
# 	showMethods("create_FAobject")
# 	Function: create_FAobject (package FAiR)
# 	restrictions="restrictions"
# 	restrictions="restrictions.1storder"
# 	restrictions="restrictions.2ndorder"
# 	restrictions="restrictions.factanal"
# 	restrictions="restrictions.general"
# 	restrictions="restrictions.orthonormal"
)
setGeneric("create_start", def = function(restrictions, ...) standardGeneric("create_start"),
		useAsDefault = TRUE, package = "FAiR"
# 	showMethods("create_start")
# 	Function: create_start (package FAiR)
# 	restrictions="restrictions"
# 	restrictions="restrictions.1storder"
# 	restrictions="restrictions.2ndorder"
# 	restrictions="restrictions.factanal"
# 	restrictions="restrictions.general"
# 	restrictions="restrictions.orthonormal"
)

## Three generics. These will be called with a user-tailored method if there are no
## specific methods available. They will not work particularly well, however.
setMethod("create_start", signature(restrictions = "restrictions"), definition = 
function(number, start, restrictions, S, model) {
	out <- apply(restrictions@Domains, 1, FUN = function(x) {
			return(runif(number, min = max(x[1], -0.75), 
					     max = min(x[2],  0.75)))
		})
	return(out)
})

setMethod("bfgs_helpS4", signature(initial = "numeric", object = "restrictions",
	done = "logical"), definition = 
function(initial, object, done = FALSE, S, lower) {
	return(NA)
})

setMethod("create_FAobject", signature(restrictions = "restrictions"), definition = 
function(restrictions, opt, manifest, call, scores, lower) {
	S <- manifest$cor
	p <- nrow(S)
	factors <- restrictions@factors
	factornames <- paste("F", 1:factors[1], sep = "")

	loadings <- array(NA_real_, dim = c(p, factors[1], 5), dimnames = c(
				rownames(S), factornames,c("PP", "PS", "RP", "RS", "FC")))

	correlations <- array(NA_real_, dim = c(factors[1], factors[1], 3), dimnames = c(
				factornames, factornames, c("PF", "RF", "PR")))


	attributes(correlations)$orthogonal <- FALSE

	trans_mats <- array(diag(factors), dim = c(factors, factors, 2), 
			dimnames = list(factornames, factornames, 
				c("primary", "reference")))


	vcov   <- matrix(NA_real_, restrictions@nvars, restrictions@nvars)
	zstats <- list()

	scores <- matrix(NA_real_, nrow = if(is.null(manifest$n.obs)) 0 else
				manifest$n.obs, ncol = factors[1])
	
	if(is.null(call$seeds)) {
		seeds <- rep(formals(Factanal)$seeds, 2)
	}
	else {
		seeds <- eval(call$seeds)
		if(length(seeds) == 1) seeds <- c(seeds, seeds)
	}
	seeds <- matrix(seeds, nrow = 1)
	colnames(seeds) <- c("unif.seed", "int.seed")

	FAobject <- new("FA",   loadings = loadings, correlations = correlations,
				trans_mats = trans_mats, uniquenesses = uniquenesses,
				restrictions = restrictions, vcov = vcov, zstats = zstats,
				scores = scores, manifest = manifest, rgenoud = list(
				extraction = opt), model = restrictions@model,
				method = restrictions@method, call = call, seeds = seeds)
				
	return(FAobject)
})

## Begin methods for restrictions.factanal
setMethod("fitS4", signature(par = "numeric", object = "restrictions.factanal"), 
definition = function(par, object, S, lower) {
	## This function is slightly modified from stats:::factanal.fit.mle, which is 
        ## Copyright 1995-2007 R Core Development Team and licensed under GPL V2+

	factors <- object@factors[1]
	sc <- 1/sqrt(par)
	Sstar <- sweep(sweep(S, 1, sc, FUN = "*"), 2, sc, FUN = "*")
	E <- eigen(Sstar, symmetric = TRUE, only.values = TRUE)
	e <- E$values[-(1:factors)]
	e <- sum(log(e) - e) - factors + nrow(S)
	return(c(fit = e))
})

setMethod("bfgs_fitS4", signature(par = "numeric", object = "restrictions.factanal",
	helper = "ANY"), definition = 
function(par, object, helper = NA, S, lower) { # same as fitS4 in this case
	## This function is slightly modified from stats:::factanal.fit.mle, which is 
        ## Copyright 1995-2007 R Core Development Team and licensed under GPL V2+

	factors <- object@factors[1]
	sc <- 1/sqrt(par)
	Sstar <- sweep(sweep(S, 1, sc, FUN = "*"), 2, sc, FUN = "*")
	E <- eigen(Sstar, symmetric = TRUE, only.values = TRUE)
	e <- E$values[-(1:factors)]
	e <- sum(log(e) - e) - factors + nrow(S)
	return(c(fit = e))
})

setMethod("gr_fitS4", signature(par = "numeric", object = "restrictions.factanal",
	helper = "ANY"),definition = 
function(par, object, helper = NA, S, lower) {
	## This function is slightly modified from stats:::factanal.fit.mle, which is 
        ## Copyright 1995-2007 R Core Development Team and licensed under GPL V2+
	q <- object@factors[1]

	Psi <- par
	sc <- 1/sqrt(Psi)
	Sstar <- sweep(sweep(S, 1, sc, FUN = "*"), 2, sc, FUN = "*")
	E <- eigen(Sstar, symmetric = TRUE)
	L <- E$vectors[, 1:q, drop = FALSE]
	load <- sweep(L, 2, sqrt(pmax(E$values[1:q] - 1, 0)), FUN = "*")
	load <- sweep(load, 1, sc, FUN = "/")
	g <- tcrossprod(load) + diag(Psi) - S
	out <- diag(g)/Psi^2
	return(out)
})

setMethod("create_start", signature(restrictions = "restrictions.factanal"), 
definition = function(restrictions, number, start, S) {
	out <- apply(restrictions@Domains, 1, FUN = function(x) {
			return(runif(number - ncol(start), min = x[1], max = x[2]))
		})
	return(rbind(t(start), out))
})

setMethod("create_FAobject", signature(restrictions = "restrictions.factanal"), 
definition = function(restrictions, opt, manifest, call, scores, lower) {

	S <- manifest$cor
	factors <- restrictions@factors[1]
	Lambda <- FAiR_get_loadings(opt$par, S, factors)
	signs <- sign(colSums(Lambda))
	signs[signs == 0] <- 1
	Lambda <- sweep(Lambda, 2, signs, FUN = "*")
	loadings <- array(Lambda, dim = c(nrow(S), factors, 5), 
			dimnames = list(rownames(S), NULL, 
					c("PP", "PS", "RP", "RS", "FC")))
	FC <- loadings[,,"FC"] <- loadings[,,"FC"]^2
	FC <- as.matrix(FC)
	sorter <- order(colSums(FC), decreasing = TRUE)
	loadings <- loadings[,sorter,,drop = FALSE]

	correlations <- array(diag(factors), dim = c(factors, factors, 3), 
			dimnames = list(NULL, NULL, c("PF", "RF", "PR")))
	attributes(correlations)$orthogonal <- factors != 1

	trans_mats <- array(diag(factors), dim = c(factors, factors, 3))

	rownames(correlations) <- colnames(correlations) <- colnames(loadings) <- 
							paste("F", 1:factors, sep = "")

	uniquenesses <- opt$par
	names(uniquenesses) <- rownames(S)

	Lambda_scaled <- sweep(Lambda, 1, sqrt(diag(manifest$cov)), FUN = "*", FALSE)
	uniquenesses_scaled <- uniquenesses * diag(manifest$cov)
	theta <- c(Lambda_scaled, uniquenesses_scaled)
	if(is.numeric(manifest$n.obs)) {
		info <- matrix(0, length(theta), length(theta))
		C_inv <- chol2inv(chol(tcrossprod(Lambda_scaled) + 
					diag(uniquenesses_scaled)))
		CiL <- C_inv %*% Lambda_scaled
		LCiL <- t(Lambda_scaled) %*% CiL
		for(y in 1:(length(Lambda))) for(z in y:length(Lambda)) {
			y_mark <- which(Lambda == Lambda[y], TRUE)
			z_mark <- which(Lambda == Lambda[z], TRUE)
			i <- y_mark[1]
			r <- y_mark[2]
			j <- z_mark[1]
			s <- z_mark[2]
			mark_row <- which(theta == Lambda_scaled[i,r])
			mark_col <- which(theta == Lambda_scaled[j,s])
			info[mark_row, mark_col] <- C_inv[i,j] * LCiL[r,s] + 
								    CiL[i,s] * CiL[j,r]
		}
		for(y in 1:length(Lambda)) for(j in 1:length(uniquenesses)) {
			y_mark <- which(Lambda == Lambda[y], TRUE)
			i <- y_mark[1]
			r <- y_mark[2]
			info[y, j + length(Lambda)] <- C_inv[i,j] * CiL[j,r]
		}
		mark <- length(Lambda)
		for(y in 1:length(uniquenesses)) for(j in y:length(uniquenesses)) {
			info[mark + y, mark + j] <- 0.5 * C_inv[y,j]^2
		}
		diag(info) <- 0.5 * diag(info)
		info <- info + t(info)
	
		g <- matrix(NA_real_, nrow = length(theta), 
				      ncol = 0.5 * ncol(Lambda) * (ncol(Lambda) - 1))
		for(y in 1:length(Lambda)) {
			y_mark <- which(Lambda == Lambda[y], TRUE)
			i <- y_mark[1]
			r <- y_mark[2]
			count <- 1
			for(u in 1:(ncol(Lambda) - 1)) for(v in (u + 1):ncol(Lambda)) {
				g[y, count] <- (r == u) * Lambda_scaled[i,v] + 
					       (r == v) * Lambda_scaled[i,u]
				g[y, count] <- g[y,count] / uniquenesses_scaled[i]
				g[i + length(Lambda), count] <- -Lambda_scaled[i,u] * 
								Lambda_scaled[i,v] / 
								uniquenesses[i]^2	    
				count <- count + 1
			}
		}
		augmented <- matrix(0,  length(theta) + ncol(g),
					length(theta) + ncol(g))
		augmented[1:length(theta),] <- cbind(info, g)
		augmented[(length(theta) + 1):nrow(augmented),1:length(theta)] <- t(g)
		vcov <- solve(augmented)[1:length(theta),1:length(theta)]
		vcov <- vcov * 2/manifest$n.obs
		zstats <- theta / sqrt(diag(vcov))
		zstats_beta <- matrix(zstats[1:length(Lambda)], nrow = nrow(Lambda), 
								ncol = ncol(Lambda))
		zstats_beta <- zstats_beta[,sorter]
		zstats_Theta2 <- zstats[(1+length(Lambda)):length(zstats)]
		zstats <- list(beta = zstats_beta, Theta2 = zstats_Theta2)
	}
	else {
		vcov <- matrix(NA, length(theta), length(theta))
		zstats <- list()
		warning("z-statistics could not be calculated because the number of",
			" observations is unknown")
	}
		
	scores <- FAiR_scores(scores, manifest$zz, Lambda, diag(factors), uniquenesses, S)
	
	if(is.null(call$seeds)) seeds <- rep(formals(Factanal)$seeds, 2)
	else {
		seeds <- eval(call$seeds)
		if(length(seeds) == 1) seeds <- c(seeds, seeds)
	}

	seeds <- matrix(seeds, nrow = 1)
	colnames(seeds) <- c("unif.seed", "int.seed")

	FA.object <- new("FA", loadings = loadings, correlations = correlations,
				trans_mats = trans_mats, uniquenesses = uniquenesses,
				restrictions = restrictions, vcov = vcov,
				zstats = zstats, scores = scores, manifest = manifest,
				rgenoud = list(extraction = opt), model = "EFA", method = "MLE", 
				call = call, seeds = seeds)
	return(FA.object)
})
## End methods for restrictions.factanal

## Begin methods for restrictions.orthonormal
setMethod("fitS4", signature(par = "numeric", object = "restrictions.orthonormal"), 
definition = function(par, object, S, lower) {

	# Make Theta2
	uniquenesses <- par[(object@beta$num_free + 1):length(par)]
	uniquenesses_min <- min(uniquenesses)
	if(uniquenesses_min < 0) { # Heywood case -> bailout
		return(c(uniquenesses_min, rep(0, 2 + length(object@criteria))))
	}
	else diag(object@Theta2$Theta2) <- uniquenesses

	# Make beta
	beta <- object@beta$beta
	beta[object@beta$free] <- par[1:object@beta$num_free]

	signs <- ifelse(colSums(beta) >= 0, 1, -1)
	if(any(signs != 1)) { # some factors do not have positive sums -> reflect
		out  <- c(1, mean(signs))
		beta <- sweep(beta, 2, signs, FUN = "*", check.margin = FALSE)
	}
	else out <- rep(1, 2)
	object@beta$beta <- beta

	return(c(out, FAiR_lexical_driver(par, object, S, lower = lower)))
})

setMethod("bfgs_helpS4", signature(initial = "numeric",  object =
	"restrictions.orthonormal", done = "logical"), definition = 
function(initial, object, done = FALSE, S, lower) {
	fits <- fitS4(initial, object, S, lower = lower)
	marker <- which(fits != 1)[1]
	return(list(fits = fits, marker = marker, done = done))
})

setMethod("bfgs_fitS4", signature(par = "numeric", object = "restrictions.orthonormal",
	helper = "ANY"), definition = 
function(par, object, helper, S, lower) {
	fits <- fitS4(par, object, S, lower = lower)
	if(helper$done) return(fits[length(fits)])
	marker <- which(fits != 1)[1]
	     if(marker < helper$marker) return(-.Machine$integer.max)
	else if(marker > helper$marker) return( .Machine$integer.max)
	else return(fits[marker])
})

setMethod("gr_fitS4", signature(par = "numeric", 
object = "restrictions.orthonormal", helper = "ANY"), definition = 
function(par, object, helper, S, lower) {
	if( (helper$marker == length(helper$fits)) && (object@method == "MLE") ) {
		# Make Theta2
		uniquenesses <- par[(object@beta$num_free + 1):length(par)]
		uniquenesses_min <- min(uniquenesses)
		if(uniquenesses_min < 0) { # Heywood case -> bailout
			return(rep(.Machine$double.eps, length(par)))
		}
		else diag(object@Theta2$Theta2) <- uniquenesses
	
		# Make beta
		beta <- object@beta$beta
		beta[object@beta$free] <- par[1:object@beta$num_free]

		signs <- ifelse(colSums(beta) >= 0, 1, -1)
		if(any(signs != 1) && !helper$done) { # negative colsums -> bailout
			return(rep(.Machine$double.eps, length(par)))
		}
		object@beta$beta <- beta

		fits <- FAiR_lexical_driver(par, object, S, lower = lower) 
		marker <- which(fits != 1)[1]
		if(marker < length(fits)) { # C not posdef -> bailout
			return(rep(.Machine$double.eps, length(par)))
		}

		gradient <- FAiR_gradient(cormat = S, do_Phi = FALSE, 
						Phi = object@Phi, 
						beta = beta, 
						Theta2 = object@Theta2$Theta2)
		gradient <- c(gradient$dF_dbeta[object@beta$free], gradient$dF_dTheta2)
	}
	else { # find numeric gradient of whatever the marginal criterion is
		gradient <- FAiR_numeric_gradient(par, object, helper, S)
	}
	return(gradient)
})


setMethod("create_start", signature(restrictions = "restrictions.orthonormal"), 
definition = function(restrictions, number, start, S) {
	factors <- restrictions@factors
	beta    <- apply(restrictions@Domains[1:sum(restrictions@beta$free),], 1, FUN =
			function(x){
				runif(number, min = -0.25, max = 0.75)
			})
	draws <- matrix(rnorm(number * nrow(S), mean = 1 - start, sd = 0.1), 
			nrow = number)
				
	beta[,1] <- ifelse(draws[,1] < 0, 0.9, 
		    ifelse(draws[,1] > 1, 0.1, sqrt(1 - draws[,1])))
	out <- cbind(beta, draws)
	return(out)
})

setMethod("create_FAobject", signature(restrictions = "restrictions.orthonormal"), 
definition = function(restrictions, opt, manifest, call, scores, lower) {
	old_restrictions <- restrictions
	S <- manifest$cor
	par <- opt$par
	factors <- restrictions@factors[1]
	fits <- opt$value
	marker <- which(fits != 1)[1]
	if(marker < length(fits)) {
		warning("constraints did not bind")
	}

	mark_start <- 1
	mark_end <- sum(restrictions@beta$free)

	Phi <- restrictions@Phi
	beta <- restrictions@beta$beta
	beta[restrictions@beta$free] <- par[mark_start:mark_end]
	Theta2 <- restrictions@Theta2$Theta2
	diag(Theta2) <- par[(mark_end + 1):length(par)]

	# construct primary pattern and factor contributions
	Pi <- beta %*% Phi
	FC <- beta * Pi
	uniquenesses <- 1 - rowSums(FC)
	names(uniquenesses) <- rownames(S)
	sorter <- order(colSums(FC), decreasing = TRUE)

	# construct reference structure and reference pattern
	Phi_inv <- chol2inv(chol(Phi))
	D <- 1/sqrt(diag(Phi_inv))
	Psi <- sweep(sweep(Phi_inv, 1, D, FUN = "*"), 2, D, FUN = "*")
	Upsilon <- sweep(beta, 2, D, FUN = "*")
	RP <- sweep(Pi, 2, D, FUN = "/")

	if( (is.numeric(manifest$n.obs)) && (restrictions@method == "MLE") ) {
		par <- sweep(beta, 1, sqrt(diag(manifest$cov)), FUN = "*", FALSE)
		par <- c(par, uniquenesses * diag(manifest$cov))
		free <- par != 0
		Hessian <- FAiR_Hessian(par, restrictions, manifest$cov, lower)
		negHessian <- -1 * Hessian[free,free]
		ev <- eigen(negHessian, TRUE, TRUE)$values
		if(ev[length(ev)] < sqrt(.Machine$double.eps)) {
			warning("Hessian indefinite, trying BHHH")
			vcov   <- FAiR_BHHH(par, old_restrictions, manifest)
		}
		else {
			vcov <- try(chol2inv(chol(negHessian)))
			if(is.matrix(vcov)) vcov <- vcov * (2/(manifest$n.obs))
			else vcov   <- FAiR_BHHH(par, old_restrictions, manifest)
		}
		zstats <- par[free] / sqrt(diag(vcov))
		zstats_beta <- beta * NA_real_
		zstats_beta[restrictions@beta$free] <- zstats[mark_start:mark_end]
		zstats_beta <- zstats_beta[,sorter]
		mark_start <- mark_end + 1
		mark_end <- length(zstats)
		zstats_Theta2 <- zstats[mark_start:mark_end]
		zstats <- list(beta = zstats_beta, Theta2 = zstats_Theta2)
	}
	else {
		vcov <- matrix(NA_real_, restrictions@nvars, restrictions@nvars)
		zstats <- list()
		warning("z-statistics could not be calculated because the number of",
			" observations is unknown")
	}

	# bundle everything up
	loadings <- array(cbind(beta, Pi, RP, Upsilon, FC),
			dim = c(nrow(S), factors, 5),
			dimnames = list(rownames(S), NULL, 
					c("PP", "PS", "RP", "RS", "FC")))
	loadings <- loadings[,sorter,]

	correlations <- array(cbind(Phi, Psi, diag(D)),
			dim = c(factors, factors, 3), 
			dimnames = list(NULL, NULL, c("PF", "RF", "PR")))
	correlations <- correlations[sorter,sorter,]
	attributes(correlations)$orthogonal <- TRUE
	rownames(correlations) <- colnames(correlations) <- colnames(loadings) <- 
							paste("F", 1:factors, sep = "")

	trans_mats <- array(diag(factors), dim = c(factors, factors, 2), 
			dimnames = list(paste("F", 1:factors), paste("F", 1:factors, 
					sep = ""), c("primary", "reference")))

	scores <- FAiR_scores(scores, manifest$zz, beta, diag(factors), uniquenesses, S)
	
	if(is.null(call$seeds)) seeds <- rep(formals(Factanal)$seeds, 2)
	else {
		seeds <- eval(call$seeds)
		if(length(seeds) == 1) seeds <- c(seeds, seeds)
	}

	seeds <- matrix(seeds, nrow = 1)
	colnames(seeds) <- c("unif.seed", "int.seed")

	FAobject <- new("FA",   loadings = loadings, correlations = correlations,
				trans_mats = trans_mats, uniquenesses = uniquenesses,
				restrictions = restrictions, vcov = vcov, zstats = zstats,
				scores = scores, manifest = manifest, 
				rgenoud = list(extraction = opt), model = "EFA", 
				method = restrictions@method, call = call, seeds = seeds)
	return(FAobject)
})
## End methods for restrictions.orthonormal

## Begin methods for restrictions.1st_order
setMethod("fitS4", signature(par = "numeric", object = "restrictions.1storder"), 
definition = function(par, object, S, lower) {
	SEFA <- object@model == "SEFA"
	# Make Theta2
	uniquenesses <- par[object@Theta2$select]
	uniquenesses_min <- min(uniquenesses)
	if(uniquenesses_min < 0) { # Heywood case -> bailout
		return(c(uniquenesses_min, rep(0, 4 + length(object@criteria))))
	}
	else diag(object@Theta2$Theta2) <- uniquenesses

	# Make Phi
	Phi <- object@Phi
	mark_start <- 1
	mark_end   <- 0.5 * object@factors[1] * (object@factors[1] - 1)
	Phi[upper.tri(Phi)] <- par[mark_start:mark_end]
	Phi <- Phi + t(Phi)
	Phi <- FAiR_nearPD(Phi, posd.tol = lower)
	object@Phi <- Phi
	check <- attributes(object@Phi)$ev
	if(check < lower) { # Phi not posdef -> bailout
		return(c(1, check, rep(0, 3 + length(object@criteria))))
	}

	# Make beta
	beta <- object@beta$beta
	beta[object@beta$free] <- par[object@beta$select]

	if(SEFA) { # enforce and check Howe (1955) conditions
		object@beta$fix_beta_args$coefs  <- beta
		object@beta$fix_beta_args$cormat <- object@Phi
		beta <- do.call(FAiR_fix_coefficients, args = object@beta$fix_beta_args)
		check <- FAiR_check_coefficients(beta, threshold = lower)
		if(check != 1) { # fail rank check -> bailout
			return(c(1, 1, check, rep(0, 2 + length(object@criteria))))
		}
	}

	signs <- ifelse(colSums(beta) >= 0, 1, -1)
	if(any(signs != 1)) { # factors not all positive -> reflect
		out  <- c(1, 1, 1, mean(signs))
		beta <- sweep(beta, MARGIN = 2, signs, FUN = "*", check.margin = FALSE)
		Phi  <- sweep(sweep(Phi, 1, signs, "*", FALSE), 2, signs, "*", FALSE)
		object@Phi <- Phi
	}
	else out <- rep(1, 4)
	object@beta$beta <- beta

	return(c(out, FAiR_lexical_driver(par, object, S, lower = lower))) 
})

setMethod("bfgs_helpS4", signature(initial = "numeric",  object =
	"restrictions.1storder", done = "logical"), definition = 
function(initial, object, done = FALSE, S, lower) {
	par <- initial
	SEFA <- object@model == "SEFA"

	# Make Theta2
	uniquenesses <- par[object@Theta2$select]
	uniquenesses_min <- min(uniquenesses)
	if(uniquenesses_min < 0) { # Heywood case -> bailout
		marker <- 1
		fits <- c(uniquenesses_min, rep(0, 4 + length(object@criteria)))
		return(list(fits = fits, marker = marker, done = done))
	}
	else diag(object@Theta2$Theta2) <- uniquenesses

	# Make Phi
	Phi <- object@Phi
	mark_start <- 1
	mark_end   <- 0.5 * object@factors[1] * (object@factors[1] - 1)
	Phi[upper.tri(Phi)] <- par[mark_start:mark_end]
	Phi <- Phi + t(Phi)
	Phi <- FAiR_nearPD(Phi, posd.tol = lower)
	object@Phi <- Phi
	check <- attributes(object@Phi)$ev
	if(check < lower) { # Phi not posdef -> bailout
		marker <- 2
		fits <- c(1, check, rep(0, 3 + length(object@criteria)))
		return(list(fits = fits, marker = marker, done = done))
	}

	# Make beta
	beta <- object@beta$beta
	beta[object@beta$free] <- par[object@beta$select]

	if(SEFA) { # enforce and check Howe (1955) conditions
		object@beta$fix_beta_args$coefs  <- beta
		object@beta$fix_beta_args$cormat <- object@Phi
		beta <- do.call(FAiR_fix_coefficients, args = object@beta$fix_beta_args)
		check <- FAiR_check_coefficients(beta, threshold = lower)
		if(check != 1) { # fail rank check -> bailout
			marker <- 3
			fits <- c(1, 1, check, rep(0, 2 + length(object@criteria)))
			return(list(fits = fits, marker = marker, done = done))
		}
	}
	object@beta$beta <- beta

	signs <- ifelse(colSums(beta) >= 0, 1, -1)
	if(any(signs != 1)) { # factors not all positive -> bailout
		marker <- 4
		fits  <- c(1, 1, 1, mean(signs), rep(0, 1 + length(object@criteria)))
		return(list(fits = fits, marker = marker, done = done))
	}
	else { # get lexical criteria
		fits <- c(rep(1, 4), FAiR_lexical_driver(par, object, S, lower = lower))
		marker <- which(fits != 1)[1]
	}
	
	if( (marker == length(fits)) && (object@model == "SEFA") ) {
		# Mark which elements of beta got squashed to zero
		squashed <- object@beta$select
		squashed[squashed] <- c(!beta) & object@beta$free
		return(list(fits = fits, marker = marker, 
			squashed = squashed, done = done))
	}
	else return(list(fits = fits, marker = marker, done = done))
})

setMethod("bfgs_fitS4", signature(par = "numeric", object = "restrictions.1storder",
	helper = "list"), definition = 
function(par, object, helper, S, lower) {
	SEFA <- object@model == "SEFA"
	if(helper$marker == length(helper$fits)) {
		if(SEFA) { # squash based on helper and pretend it's a CFA
			object@model <- "CFA"
			small <- par[helper$squashed]
			par[helper$squashed] <- 0
		}
		fits <- fitS4(par, object, S, lower = lower) # get lexical criteria
		if(helper$done) return(fits[length(fits)])
		marker <- which(fits != 1)[1]
		if(marker < length(helper$fits)) return(-.Machine$double.xmax)

		if(SEFA) out <- fits[marker] - crossprod(small) # penalize if small != 0
		else     out <- fits[marker]
		return(out)
	}
	else {
		fits <- fitS4(par, object, S, lower = lower) # get lexical criteria
		marker <- which(fits != 1)[1]
		if(marker > helper$marker) return(.Machine$double.xmax)
		return(fits[marker])
	}
})

setMethod("gr_fitS4", signature(par = "numeric", object = "restrictions.1storder", 
	helper = "list"), definition = 
function(par, object, helper, S, lower) {
	SEFA <- object@model == "SEFA"
	if( (helper$marker == length(helper$fits)) && (object@method == "MLE") ) {
		# Make Theta2
		uniquenesses <- par[object@Theta2$select]
		uniquenesses_min <- min(uniquenesses)
		if(uniquenesses_min < 0) { # Heywood case -> bailout
			return(rep(.Machine$double.eps, length(par)))
		}
		else diag(object@Theta2$Theta2) <- uniquenesses

		# Make Phi
		Phi <- object@Phi
		mark_start <- 1
		mark_end   <- 0.5 * object@factors[1] * (object@factors[1] - 1)
		Phi[upper.tri(Phi)] <- par[mark_start:mark_end]
		Phi <- Phi + t(Phi)
		Phi <- FAiR_nearPD(Phi, posd.tol = lower)
		object@Phi <- Phi
		check <- attributes(Phi)$ev
		if(check < lower) { # Phi not posdef -> bailout
			return(rep(.Machine$double.eps, length(par)))
		}

		# Make beta
		beta <- object@beta$beta
	
		if(SEFA) { # squash based on helper and pretend it's a CFA
			small <- par[helper$squashed]
			par[helper$squashed] <- 0
		}

		beta[object@beta$free] <- par[object@beta$select]

		signs <- ifelse(colSums(beta) >= 0, 1, -1)
		if(any(signs != 1) && !helper$done) { # factors not positive -> bailout
			return(rep(.Machine$double.eps, length(par)))
		}
		object@beta$beta <- beta

		fits <- FAiR_lexical_driver(par, object, S, lower = lower) 
		marker <- which(fits != 1)[1]
		if(marker < length(fits)) { # failed a constraint -> bailout 
			return(rep(.Machine$double.eps, length(par)))
		}

		gradient <- FAiR_gradient(cormat = S, do_Phi = TRUE, 
						Phi = Phi, 
						beta = beta, 
						Theta2 = object@Theta2$Theta2)
		gradient <- c(  gradient$dF_dPhi[upper.tri(Phi)] * 2,
				gradient$dF_dbeta[object@beta$free], 
				gradient$dF_dTheta2)
		if(SEFA) { # quadratic loss around zero for small coefficients if !done
			gradient[helper$squashed] <- -2 * small * (!helper$done)
		}
	}
	else { # find numeric gradient of whatever the marginal criterion is
		gradient <- FAiR_numeric_gradient(par, object, helper, S)
	}
	return(gradient)
})

setMethod("create_start", signature(restrictions = "restrictions.1storder"), 
definition = function(restrictions, number, start, S) {

	if(restrictions@model == "CFA") {
		out <- apply(restrictions@Domains, 1, FUN = function(x) {
				return(runif(number, min = max(x[1], -0.75), 
						     max = min(x[2],  0.75)))
			})
		return(out)		
	}
	factors <- restrictions@factors[1]
	Lambda  <- FAiR_get_loadings(1 - c(start), S, factors)
	Lambda_swept <- sweep(Lambda, 1, sqrt(rowSums(Lambda^2)), FUN = "/", FALSE)
	combos <- combn(1:nrow(Lambda), factors)
	out  <- matrix(NA_real_, nrow = ncol(combos), ncol = restrictions@nvars)
	fits <- matrix(NA_real_, nrow = ncol(combos), ncol = 5 + 
			length(restrictions@criteria))
	for(i in 1:ncol(combos)) {
		Tmat <- t(Lambda_swept[combos[,i],])
		beta  <- try(Lambda %*% t(solve(Tmat)))
		if(!is.matrix(beta)) next
		signs <- sign(colSums(beta))
		if(any(signs != 1)) {
			beta <- sweep(beta, 2, signs, FUN = "*", FALSE)
			Phi_temp <- crossprod(Tmat)
			Phi_temp <- sweep(sweep(Phi_temp, 1, signs, FUN = "*", FALSE),
							  2, signs, FUN = "*", FALSE)
		}
		else Phi_temp <- crossprod(Tmat)
		Phi_temp <- FAiR_nearPD(Phi_temp)
		FC <- beta * (beta %*% Phi_temp)
		sorter <- order(colSums(FC), decreasing = TRUE)
		beta <- beta[,sorter]
		Phi_temp <- Phi_temp[sorter, sorter]
		uniquenesses <- rnorm(nrow(start), 1 - start, sd = .1)
		out[i,] <- c(Phi_temp[upper.tri(Phi_temp)], 
				c(beta[restrictions@beta$free]), uniquenesses)
		out[i,] <-      ifelse(out[i,] < restrictions@Domains[,1], 
						 restrictions@Domains[,1],
				ifelse(out[i,] > restrictions@Domains[,2], 
						 restrictions@Domains[,2], out[i,]))
		fits[i,] <- fitS4(out[i,], restrictions, S, 
				    lower = sqrt(.Machine$double.eps))
	}
	out <- cbind(out, fits)
	out <- out[complete.cases(out),]
	out <- out[do.call(order, as.data.frame(fits[complete.cases(fits),])),]
	return(out)
})

setMethod("create_FAobject", signature(restrictions = "restrictions.1storder"),
definition = function(restrictions, opt, manifest, call, scores, lower) {
	old_restrictions <- restrictions
	SEFA <- restrictions@model == "SEFA"
	S <- manifest$cor
	old_par <- par <- opt$par
	factors <- restrictions@factors[1]
	fits <- opt$value
	marker <- which(fits != 1)[1]
	if(marker < length(fits)) {
		warning("constraints did not bind")
	}

	# Make Theta2
	uniquenesses <- par[restrictions@Theta2$select]

	# Make Phi
	Phi <- restrictions@Phi
	mark_start <- 1
	mark_end   <- 0.5 * factors * (factors - 1)
	Phi[upper.tri(Phi)] <- par[mark_start:mark_end]
	Phi <- Phi + t(Phi)
	Phi <- FAiR_nearPD(Phi, posd.tol = lower)

	# Make beta
	beta <- restrictions@beta$beta
	beta[restrictions@beta$free] <- par[restrictions@beta$select]

	if(SEFA) { # enforce and check Howe (1955) conditions
		restrictions@beta$fix_beta_args$coefs  <- beta
		restrictions@beta$fix_beta_args$cormat <- Phi
		beta <- do.call(FAiR_fix_coefficients, 
				args = restrictions@beta$fix_beta_args)
	}

	# construct primary pattern and factor contributions
	Pi <- beta %*% Phi
	FC <- beta * Pi
	uniquenesses <- 1 - rowSums(FC)
	names(uniquenesses) <- rownames(S)
	sorter <- order(colSums(FC), decreasing = TRUE)

	# construct reference structure and reference pattern
	Phi_inv <- chol2inv(chol(Phi))
	D <- 1/sqrt(diag(Phi_inv))
	Psi <- sweep(sweep(Phi_inv, 1, D, FUN = "*"), 2, D, FUN = "*")
	Upsilon <- sweep(beta, 2, D, FUN = "*")
	RP <- sweep(Pi, 2, D, FUN = "/")

	loadings <- array(cbind(beta, Pi, RP, Upsilon, FC),
			dim = c(nrow(S), factors, 5),
			dimnames = list(rownames(S), NULL, 
					c("PP", "PS", "RP", "RS", "FC")))

	correlations <- array(cbind(Phi, Psi, diag(D)), dim = c(factors, factors, 3),
			dimnames = list(NULL, NULL, c("PF", "RF", "PR")))

	trans_mats <- array(diag(factors), dim = c(factors, factors, 2), 
			dimnames = list(paste("F", 1:factors), paste("F", 1:factors, 
					sep = ""), c("primary", "reference")))

	if(is.numeric(manifest$n.obs)) {
		par <- sweep(loadings[,,1], 1, sqrt(diag(manifest$cov)), FUN = "*", FALSE)
		par <- c(Phi[upper.tri(Phi)], par[restrictions@beta$free], 
			uniquenesses * diag(manifest$cov))
		free <- par != 0
		Hessian <- FAiR_Hessian(par, old_restrictions, manifest$cov, lower)
		negHessian <- -1 * Hessian[free,free]
		ev <- eigen(negHessian, TRUE, TRUE)$values
		if(ev[length(ev)] < 0) {
			warning("Hessian indefinite, trying BHHH")
			vcov   <- FAiR_BHHH(par, old_restrictions, manifest)
		}
		else {
			vcov <- try(chol2inv(chol(negHessian)))
			if(is.matrix(vcov)) vcov <- vcov * (2/(manifest$n.obs))
			else vcov   <- FAiR_BHHH(par, old_restrictions, manifest)
		}
		zstats  <- par[free] / sqrt(diag(vcov))
		zstats_Phi  <- 0 * Phi
		mark_start <- 1
		mark_end <- 0.5 * factors * (factors - 1)
		zstats_Phi[upper.tri(zstats_Phi)] <- zstats[mark_start:mark_end]
		zstats_Phi <- zstats_Phi + t(zstats_Phi)
		diag(zstats_Phi) <- NA_real_
		zstats_Phi <- zstats_Phi[sorter,sorter]
		zstats_beta <- beta * NA_real_
		mark_start <- mark_end + 1
		mark_end <- length(zstats) - nrow(S)
		fill <- restrictions@beta$free & c(beta != 0)
		zstats_beta[fill] <- zstats[mark_start:mark_end]
		zstats_beta <- zstats_beta[,sorter]
		mark_start <- mark_end + 1
		mark_end <- length(zstats)
		zstats_Theta2 <- zstats[mark_start:mark_end]
		zstats <- list(Phi = zstats_Phi, beta = zstats_beta, 
				Theta2 = zstats_Theta2)
	}
	else {
		warning("z-statistics could not be calculated")
		vcov <- matrix(NA_real_,nrow = restrictions@nvars, 
					ncol = restrictions@nvars)
		zstats <- list()
	}

	loadings <- loadings[,sorter,]
	correlations <- correlations[sorter,sorter,]
	attributes(correlations)$orthogonal <- FALSE
	parnames <- paste("F", 1:factors, sep = "")
	rownames(correlations) <- colnames(correlations) <- colnames(loadings) <- parnames

	scores <- FAiR_scores(scores, manifest$zz, beta, Phi, uniquenesses, S)
	
	if(is.null(call$seeds)) seeds <- rep(formals(Factanal)$seeds, 2)
	else {
		seeds <- eval(call$seeds)
		if(length(seeds) == 1) seeds <- c(seeds, seeds)
	}
	seeds <- matrix(seeds, nrow = 1)
	colnames(seeds) <- c("unif.seed", "int.seed")

	if(SEFA) { # recalculate degrees of freedom
		dof <- 0.5 * ncol(S) * (ncol(S) + 1) - restrictions@nvars
		dof <- dof + sum(restrictions@Domains[,1] == restrictions@Domains[,2])
		dof <- dof + sum(!beta)
		dof <- dof + sum(restrictions@beta$beta != 0, na.rm = TRUE)
		restrictions@dof <- as.integer(dof)
	}
	FAobject <- new("FA",   loadings = loadings, correlations = correlations,
				trans_mats = trans_mats, uniquenesses = uniquenesses,
				restrictions = restrictions, vcov = vcov, zstats = zstats,
				scores = scores, manifest = manifest, 
				rgenoud = list(extraction = opt), model = restrictions@model, 
				method = restrictions@method, call = call, seeds = seeds)
	return(FAobject)
})
## end methods for restrictions.1storder


## begin methods for restrictions.general
setMethod("fitS4", signature(par = "numeric", object = "restrictions.general"), 
definition = function(par, object, S, lower) {
	SEFA <- object@model == "SEFA"

	# Make Theta2
	uniquenesses <- par[object@Theta2$select]
	uniquenesses_min <- min(uniquenesses)
	if(uniquenesses_min < 0) { # Heywood case -> bailout
		return(c(uniquenesses_min, rep(0, 4 + length(object@criteria))))
	}
	else diag(object@Theta2$Theta2) <- uniquenesses

	# Make Delta and Phi
	object@Delta$Delta[object@Delta$free] <- par[1:object@Delta$num_free]
	Phi   <- tcrossprod(object@Delta$Delta)
	diag(Phi) <- 1
	check <- eigen(Phi, TRUE, TRUE)$values[ncol(Phi)]
	if(check < lower) { # Phi not posdef -> bailout
		return(c(1, check, rep(0, 3 + length(object@criteria))))
	}
	object@Phi <- Phi

	# Make beta
	beta <- object@beta$beta
	beta[object@beta$free] <- par[object@beta$select]

	if(SEFA) { # enforce and check Howe (1955) conditions
		object@beta$fix_beta_args$coefs  <- beta
		object@beta$fix_beta_args$cormat <- object@Phi
		beta <- do.call(FAiR_fix_coefficients, args = object@beta$fix_beta_args)
		check <- FAiR_check_coefficients(beta, threshold = lower)
		if(check != 1) { # fail rank check -> flag this
			out <- c(1, 1, check)
		}
		else out <- rep(1,3)
	}
	out <- rep(1,3)

	signs <- ifelse(colSums(beta) >= 0, 1, -1)
	Delta_sum <- sum(object@Delta$Delta)
	if(any(signs != 1)) { # factors not all positive -> reflect
		out  <- c(out, (sum(signs) + (Delta_sum >= 0)) / (length(signs) + 1))
		beta <- sweep(beta, MARGIN = 2, signs, FUN = "*", check.margin = FALSE)
		Phi  <- sweep(sweep(Phi, 1, signs, "*", FALSE), 2, signs, "*", FALSE)
		object@Phi <- Phi
	}
	else if(Delta_sum < 0) {
		object@Delta$Delta <- -1 * object@Delta$Delta
		out <- c(out, length(signs) / (length(signs) + 1))
	}
	else out <- c(out, 1)
	object@beta$beta <- beta

	return(c(out, FAiR_lexical_driver(par, object, S, Delta, lower = lower)))
})

setMethod("bfgs_helpS4", signature(initial = "numeric",  object = "restrictions.general",
	done = "logical"), definition = 
function(initial, object, done = FALSE, S, lower) {
	par <- initial
	SEFA <- object@model == "SEFA"

	# Make Theta2
	uniquenesses <- par[object@Theta2$select]
	uniquenesses_min <- min(uniquenesses)
	if(uniquenesses_min < 0) { # Heywood case -> bailout
		marker <- 1
		fits <- c(uniquenesses_min, rep(0, 4 + length(object@criteria)))
		return(list(fits = fits, marker = marker, done = done))
	}
	else diag(object@Theta2$Theta2) <- uniquenesses

	# Make Delta and Phi
	object@Delta$Delta[object@Delta$free] <- par[1:object@Delta$num_free]
	Phi   <- tcrossprod(object@Delta$Delta)
	diag(Phi) <- 1
	check <- eigen(Phi, TRUE, TRUE)$values[ncol(Phi)]
	if(check < lower) { # Phi not posdef -> bailout
		marker <- 2
		fits <- c(1, check, rep(0, 3 + length(object@criteria)))
		return(list(fits = fits, marker = marker, done = done))
	}
	object@Phi <- Phi

	# Make beta
	beta <- object@beta$beta
	beta[object@beta$free] <- par[object@beta$select]

	if(SEFA) { # enforce and check Howe (1955) conditions
		object@beta$fix_beta_args$coefs  <- beta
		object@beta$fix_beta_args$cormat <- object@Phi
		beta <- do.call(FAiR_fix_coefficients, args = object@beta$fix_beta_args)
		check <- FAiR_check_coefficients(beta, threshold = lower)
		if(check != 1) { # fail rank check -> bailout
			marker <- 3
			fits <- c(1, 1, check, rep(0, 2 + length(object@criteria)))
			return(list(fits = fits, marker = marker, done = done))
		}
	}
	object@beta$beta <- beta

	Delta_sum <- sum(object@Delta$Delta)
	signs <- ifelse(colSums(beta) >= 0, 1, -1)
	if(any(signs != 1)) { # factors not all positive -> bailout
		marker <- 4
		fits  <- c(1, 1, 1, (sum(signs) + (Delta_sum >= 0)) / 
			(length(signs) + 1), rep(0, 1 + length(object@criteria)))
		return(list(fits = fits, marker = marker, done = done))
	}
	else if(Delta_sum < 0) {
		marker <- 4
		fits <- c(1, 1, 1, length(signs) / (length(signs) + 1), 
				rep(0, 1 + length(object@criteria)))
		return(list(fits = fits, marker = marker, done = done))
	}
	else { # get lexical criteria
		fits <- c(rep(1, 4), FAiR_lexical_driver(par, object, S, 
							Delta, lower = lower))
		marker <- which(fits != 1)[1]
	}
	
	if( (marker == length(fits)) && SEFA ) {
		# Mark which elements of beta got squashed to zero
		squashed <- object@beta$select
		squashed[squashed] <- c(!beta) & object@beta$free
		return(list(fits = fits, marker = marker, 
			squashed = squashed, done = done))
	}
	else return(list(fits = fits, marker = marker, done = done))
})

setMethod("bfgs_fitS4", signature(par = "numeric", object = "restrictions.general", 
	helper = "list"), definition = 
function(par, object, helper, S, lower) {
	SEFA <- object@model == "SEFA"
	if(helper$marker == length(helper$fits)) {
		if(SEFA) { # squash based on helper and pretend it's a CFA
			object@model <- "CFA"
			small <- par[helper$squashed]
			par[helper$squashed] <- 0
		}
		fits <- fitS4(par, object, S, lower = lower) # get lexical criteria
		if(helper$done) return(fits[length(fits)])
		marker <- which(fits != 1)[1]
		if(marker < length(helper$fits)) return(-.Machine$integer.max)

		if(SEFA) out <- fits[marker] - crossprod(small) # penalize if small != 0
		else     out <- fits[marker]
		return(out)
	}
	else {
		fits <- fitS4(par, object, S, lower = lower) # get lexical criteria
		marker <- which(fits != 1)[1]
		if(marker > helper$marker) return(.Machine$integer.max)
		return(fits[marker])
	}
})

setMethod("gr_fitS4", signature(par = "numeric", object = "restrictions.general", 
	helper = "list"), definition = 
function(par, object, helper, S, lower) {
	SEFA <- object@model == "SEFA"
	if( (helper$marker == length(helper$fits)) && (object@method == "MLE") ) {
		# Make Theta2
		uniquenesses <- par[object@Theta2$select]
		uniquenesses_min <- min(uniquenesses)
		if(uniquenesses_min < 0) { # Heywood case -> bailout
			return(rep(.Machine$double.eps, length(par)))
		}
		else diag(object@Theta2$Theta2) <- uniquenesses

		# Make Delta and Phi
		object@Delta$Delta[object@Delta$free] <- par[1:object@Delta$num_free]
		if(sum(object@Delta$Delta) < 0) { # Delta negative-sum -> bailout
			return(rep(.Machine$double.eps, length(par)))
		}
		Phi   <- tcrossprod(object@Delta$Delta)
		diag(Phi) <- 1
		check <- eigen(Phi, TRUE, TRUE)$values[ncol(Phi)]
		if(check < lower) { # Phi not posdef -> bailout
			return(rep(.Machine$double.eps, length(par)))
		}
		object@Phi <- Phi

		# Make beta
		beta <- object@beta$beta
	
		if(SEFA) { # squash based on helper and pretend it's a CFA
			small <- par[helper$squashed]
			par[helper$squashed] <- 0
		}

		beta[object@beta$free] <- par[object@beta$select]

		signs <- ifelse(colSums(beta) >= 0, 1, -1)
		if(any(signs != 1) && !helper$done) { # factors not positive -> bailout
			return(rep(.Machine$double.eps, length(par)))
		}
		object@beta$beta <- beta

		fits <- FAiR_lexical_driver(par, object, S, Delta, lower = lower)
		marker <- which(fits != 1)[1]
		if(marker < length(fits)) { # failed a constraint -> bailout 
			return(rep(.Machine$double.eps, length(par)))
		}

		gradient <- FAiR_gradient(cormat = S, do_Phi = TRUE, 
						Phi = Phi, 
						beta = beta, 
						Theta2 = object@Theta2$Theta2)
		dF_dPhi <- gradient$dF_dPhi
		diag(dF_dPhi) <- 0
		dF_dDelta <- rep(0, length(object@Delta$Delta))
		for(i in 1:length(dF_dDelta)) dF_dDelta[i] <- sum(2 * dF_dPhi[i,] * 
								object@Delta$Delta)
		dF_dDelta <- dF_dDelta[object@Delta$free]
		gradient <- c(dF_dDelta, gradient$dF_dbeta[object@beta$free], 
				gradient$dF_dTheta2)

		if(SEFA) { # quadratic loss around zero for small coefficients if !done
			gradient[helper$squashed] <- -2 * small * (!helper$done)
		}
	}
	else { # find numeric gradient of whatever the marginal criterion is
		gradient <- FAiR_numeric_gradient(par, object, helper, S)
	}
	return(gradient)
})

setMethod("create_start", signature(restrictions = "restrictions.general"), 
	definition = function(restrictions, number, start, S) {
	factors <- restrictions@factors[1]
	if(restrictions@model == "CFA") {
		out <- apply(restrictions@Domains, 1, FUN = function(x) {
				return(runif(number, min = max(x[1], -0.75), 
						     max = min(x[2],  0.75)))
			})
		return(out)		
	}
	Lambda  <- FAiR_get_loadings(1 - c(start), S, factors)
	signs <- sign(colSums(Lambda))
	if(any(signs != 1)) Lambda <- sweep(Lambda, 2, signs, FUN = "*", FALSE)
	Lambda <- Lambda %*% t(solve(FAiR_Landahl(factors)))
	out  <- matrix(NA_real_, nrow = number, ncol = restrictions@nvars)
	while(number) {
		Delta <- as.matrix(runif(factors, min = restrictions@Domains[1:factors,1],
						  max = restrictions@Domains[1:factors,2]))
		if(sum(Delta) < 0) Delta <- Delta * -1
		Phi <- tcrossprod(Delta)
		diag(Phi) <- 1
		Tmat <- try(chol(Phi), silent = TRUE)
		if(!is.matrix(Tmat)) next
		beta <- Lambda %*% t(chol2inv(Tmat))
		signs <- sign(colSums(beta))
		if(any(signs != 1)) beta <- sweep(beta, 2, signs, FUN = "*", FALSE)
		FC <- beta * (beta %*% Phi)
		sorter <- order(colSums(FC), decreasing = TRUE)
		beta <- beta[,sorter]
		Delta <- Delta[sorter,]
		uniquenesses <- rnorm(nrow(start), 1 - start, sd = .1)
		out[number,] <- c(Delta, beta[restrictions@beta$free], uniquenesses)
		out[number,] <- ifelse(out[number,] < restrictions@Domains[,1], 
						      restrictions@Domains[,1],
				ifelse(out[number,] > restrictions@Domains[,2], 
						      restrictions@Domains[,2], 
						      out[number,]))
		number <- number - 1
	}
	return(out)
})

setMethod("create_FAobject", signature(restrictions = "restrictions.general"), 
definition = function(restrictions, opt, manifest, call, scores, lower) {
	old_restrictions <- restrictions
	SEFA <- restrictions@model == "SEFA"
	S <- manifest$cor
	par <- opt$par
	factors <- restrictions@factors[1]
	fits <- opt$value
	marker <- which(fits != 1)[1]
	if(marker < length(fits)) {
		warning("constraints did not bind")
	}

	# Make Theta2
	uniquenesses <- par[restrictions@Theta2$select]
	Theta2 <- diag(uniquenesses)

	# Make Delta and Phi
	Delta <- restrictions@Delta$Delta
	Delta[restrictions@Delta$free] <- par[1:restrictions@Delta$num_free]
	Phi   <- tcrossprod(Delta)
	diag(Phi) <- 1

	# Make beta
	beta <- restrictions@beta$beta
	beta[restrictions@beta$free] <- par[restrictions@beta$select]

	if(SEFA) { # enforce and check Howe (1955) conditions
		restrictions@beta$fix_beta_args$coefs  <- beta
		restrictions@beta$fix_beta_args$cormat <- Phi
		beta <- do.call(FAiR_fix_coefficients, 
				args = restrictions@beta$fix_beta_args)
	}

	uniquenesses_2nd <- c(1 - Delta^2)
	loadings_2nd <- Delta

	# construct primary pattern and factor contributions
	Pi <- beta %*% Phi
	FC <- beta * Pi
	uniquenesses <- 1 - rowSums(FC)
	names(uniquenesses) <- rownames(S)
	sorter <- order(colSums(FC), decreasing = TRUE)

	# construct reference structure and reference pattern
	Phi_inv <- chol2inv(chol(Phi))
	D <- 1/sqrt(diag(Phi_inv))
	Psi <- sweep(sweep(Phi_inv, 1, D, FUN = "*"), 2, D, FUN = "*")
	Upsilon <- sweep(beta, 2, D, FUN = "*")
	RP <- sweep(Pi, 2, D, FUN = "/")

	loadings <- array(cbind(beta, Pi, RP, Upsilon, FC),
			dim = c(nrow(S), factors, 5),
			dimnames = list(rownames(S), NULL, 
					c("PP", "PS", "RP", "RS", "FC")))

	correlations <- array(cbind(Phi, Psi, diag(D)),
			dim = c(factors, factors, 3), 
			dimnames = list(NULL, NULL, c("PF", "RF", "PR")))

	trans_mats <- array(diag(factors), dim = c(factors, factors, 2), 
			dimnames = list(paste("F", 1:factors), paste("F", 1:factors, 
					sep = ""), c("primary", "reference")))

	if(is.numeric(manifest$n.obs)) {
		par <- sweep(loadings[,,1], 1, sqrt(diag(manifest$cov)), FUN = "*", FALSE)
		par <- c(Delta[restrictions@Delta$free], par[restrictions@beta$free], 
			uniquenesses * diag(manifest$cov))
		free <- par != 0
		Hessian <- FAiR_Hessian(par, old_restrictions, manifest$cov, lower)
		negHessian <- -1 * Hessian[free,free]
		ev <- eigen(negHessian, TRUE, TRUE)$values
		if(ev[length(ev)] < 0) {
			warning("Hessian indefinite, trying BHHH estimator")
			vcov   <- FAiR_BHHH(par, old_restrictions, manifest)
		}
		else {
			vcov <- try(chol2inv(chol(negHessian)))
			if(is.matrix(vcov)) vcov <- vcov * (2/(manifest$n.obs))
			else vcov <- FAiR_BHHH(par, old_restrictions, manifest)
		}
		zstats  <- par[free] / sqrt(diag(vcov))
		mark_start <- 1
		mark_end <- factors
		zstats_Delta  <- zstats[mark_start:mark_end]
		zstats_Delta <- zstats_Delta[sorter]
		mark_start <- mark_end + 1
		mark_end <- length(zstats) - nrow(S)
		zstats_beta <- beta * NA_real_
		fill <- restrictions@beta$free & c(beta != 0)
		zstats_beta[fill] <- zstats[mark_start:mark_end]
		zstats_beta <- zstats_beta[,sorter]
		mark_start <- mark_end + 1
		mark_end <- length(zstats)
		zstats_Theta2 <- zstats[mark_start:mark_end]
		zstats <- list(Delta = zstats_Delta, beta = zstats_beta, 
				Theta2 = zstats_Theta2)
	}
	else {
		warning("zstatistics could not be calculated")
		vcov <- matrix(NA_real_,nrow = restrictions@nvars, 
					ncol = restrictions@nvars)
		zstats <- list()
	}
	loadings <- loadings[,sorter,]
	correlations <- correlations[sorter,sorter,]
	attributes(correlations)$orthogonal <- FALSE
	rownames(correlations) <- colnames(correlations) <- colnames(loadings) <- 
							paste("F", 1:factors, sep = "")

	scores <- FAiR_scores(scores, manifest$zz, beta, Phi, uniquenesses, S)
	
	if(is.null(call$seeds)) seeds <- rep(formals(Factanal)$seeds, 2)
	else {
		seeds <- eval(call$seeds)
		if(length(seeds) == 1) seeds <- c(seeds, seeds)
	}

	seeds <- matrix(seeds, nrow = 1)
	colnames(seeds) <- c("unif.seed", "int.seed")

	loadings_2nd <- loadings_2nd[sorter,, drop = FALSE]
	if(SEFA) { # recalculate degrees of freedom
		dof <- 0.5 * ncol(S) * (ncol(S) + 1) - restrictions@nvars
		dof <- dof + sum(restrictions@Domains[,1] == restrictions@Domains[,2])
		dof <- dof + sum(!beta)
		dof <- dof + sum(restrictions@beta$beta   != 0, na.rm = TRUE)
		dof <- dof + sum(restrictions@Delta$Delta != 0, na.rm = TRUE)
		restrictions@dof <- as.integer(dof)
	}

	FAobject <- new("FA.general", loadings = loadings, 
			loadings_2nd = loadings_2nd, correlations = correlations,
			trans_mats = trans_mats, uniquenesses = uniquenesses,
			uniquenesses_2nd = uniquenesses_2nd,
			restrictions = restrictions, vcov = vcov, zstats = zstats,
			scores = scores, manifest = manifest, rgenoud = list(extraction = opt), 
			model = restrictions@model, method = restrictions@method, 
			call = call, seeds = seeds)
	return(FAobject)
})
## end methods for restrictions.general

## begin methods for restrictions.2ndorder
setMethod("fitS4", signature(par = "numeric", object = "restrictions.2ndorder"), 
definition = function(par, object, S, lower) {
	SEFA <- object@model == "SEFA"

	# Make Theta2
	uniquenesses <- par[object@Theta2$select]
	uniquenesses_min <- min(uniquenesses)
	if(uniquenesses_min < 0) { # Heywood case -> bailout
		return(c(uniquenesses_min, rep(0, 8 + length(object@criteria))))
	}
	else diag(object@Theta2$Theta2) <- uniquenesses

	# Make Xi
	mark_start <- 1
	mark_end   <- 0.5 * object@factors[2] * (object@factors[2] - 1)
	Xi <- object@Xi
	Xi[upper.tri(Xi)] <- par[mark_start:mark_end]
	Xi <- Xi + t(Xi)
	Xi <- FAiR_nearPD(Xi, posd.tol = lower)
	object@Xi <- Xi
	if(attributes(Xi)$ev < lower) {
		return(c(1, attributes(object@Xi)$ev, 
			rep(0, 7 + length(object@criteria))))
	}

	# Make Delta
	Delta <- object@Delta$Delta
	Delta[object@Delta$free] <- par[object@Delta$select]
	if(SEFA) { # enforce and check Howe (1955) conditions on Delta
		object@Delta$fix_Delta_args$coefs <- Delta
		object@Delta$fix_Delta_args$cormat <- Xi
		Delta <- do.call(FAiR_fix_coefficients, args =object@Delta$fix_Delta_args)
		check <- FAiR_check_coefficients(Delta, threshold = lower)
		if(check != 1) { # fail rank check -> bailout
			return(c(1, 1, check, rep(0, 6 + length(object@criteria))))
		}
	}

	signs <- ifelse(colSums(Delta) >= 0, 1, -1)
	if(any(signs != 1)) { # 2nd order factors not all positive -> reflect
		out  <- c(1, 1, 1, mean(signs))
		Delta <- sweep(Delta, MARGIN = 2, signs, FUN = "*", check.margin = FALSE)
		Xi    <- sweep(sweep(Xi, 1, signs, "*", FALSE), 2, signs, "*", FALSE)
	}
	else out <- rep(1, 4)

	# Make Phi
	Phi <- crossprod(chol(Xi) %*% t(Delta))
	check <- max(diag(Phi))
	if(check > 1) { # Heywood case -> bailout
		return(c(out, -check, rep(0, 4 + length(object@criteria))))
	}
	diag(Phi) <-1
	object@Phi <- Phi
	out <- c(out, 1)

	check <- eigen(Phi, TRUE, TRUE)$values[ncol(Phi)]
	if(check < lower) { # Phi not posdef -> bailout
		return(c(out, check, rep(0, 3 + length(object@criteria))))
	}
	out <- c(out, 1)

	# Make beta
	beta <- object@beta$beta
	beta[object@beta$free] <- par[object@beta$select]

	if(SEFA) { # enforce and check Howe (1955) conditions on beta
		object@beta$fix_beta_args$coefs  <- beta
		object@beta$fix_beta_args$cormat <- object@Phi
		beta <- do.call(FAiR_fix_coefficients, args = object@beta$fix_beta_args)
		check <- FAiR_check_coefficients(beta, threshold = lower)
		if(check != 1) { # fail rank check -> bailout
			return(c(out, check, rep(0, 2 + length(object@criteria))))
		}
	}
	out <- c(out, 1)

	signs <- ifelse(colSums(beta) >= 0, 1, -1)
	if(any(signs != 1)) { # factors not all positive -> reflect
		out  <- c(out, mean(signs))
		beta <- sweep(beta, MARGIN = 2, signs, FUN = "*", check.margin = FALSE)
		Phi  <- sweep(sweep(Phi, 1, signs, "*", FALSE), 2, signs, "*", FALSE)
		object@Phi <- Phi
	}
	else out <- c(out, 1)
	object@beta$beta <- beta

	return(c(out, FAiR_lexical_driver(par, object, S, Delta, Xi, lower = lower))) 
})

setMethod("bfgs_helpS4", signature(initial = "numeric",  object =
	"restrictions.2ndorder", done = "logical"), definition = 
function(initial, object, done = FALSE, S, lower) {
	par <- initial
	SEFA <- object@model == "SEFA"

	# Make Theta2
	uniquenesses <- par[object@Theta2$select]
	uniquenesses_min <- min(uniquenesses)
	if(uniquenesses_min < 0) { # Heywood case -> bailout
		marker <- 1
		fits <- c(uniquenesses_min, rep(0, 8 + length(object@criteria)))
		return(list(fits = fits, marker = marker, done = done))
	}
	else diag(object@Theta2$Theta2) <- uniquenesses

	# Make Xi
	mark_start <- 1
	mark_end   <- 0.5 * object@factors[2] * (object@factors[2] - 1)
	Xi <- object@Xi
	Xi[upper.tri(Xi)] <- par[mark_start:mark_end]
	Xi <- Xi + t(Xi)
	Xi <- FAiR_nearPD(Xi, posd.tol = lower)
	object@Xi <- Xi
	if(attributes(Xi)$ev < lower) {
		marker <- 2
		fits <- c(1, attributes(object@Xi)$ev, 
				rep(0, 7 + length(object@criteria)))
		return(list(fits = fits, marker = marker, done = done))
	}

	# Make Delta
	Delta <- object@Delta$Delta
	Delta[object@Delta$free] <- par[object@Delta$select]
	if(SEFA) { # enforce and check Howe (1955) conditions on Delta
		object@Delta$fix_Delta_args$coefs  <- Delta
		object@Delta$fix_Delta_args$cormat <- Xi
		Delta <- do.call(FAiR_fix_coefficients, args =object@Delta$fix_Delta_args)
		check <- FAiR_check_coefficients(Delta, threshold = lower)
		if(check != 1) { # fail rank check -> bailout
			marker <- 3
			fits <- c(1, 1, check, rep(0, 6 + length(object@criteria)))
			return(list(fits = fits, marker = marker, done = done))
		}
	}

	signs <- ifelse(colSums(Delta) >= 0, 1, -1)
	if(any(signs != 1) && !done) { # 2nd order factors not all positive -> bailout
		marker <- 4
		fits   <- c(rep(1,3), mean(signs), rep(0, 5 + length(object@criteria)))
		return(list(fits = fits, marker = marker, done = done))
	}

	# Make Phi
	Phi <- crossprod(chol(Xi) %*% t(Delta))
	check <- max(diag(Phi))
	if(check > 1) { # Heywood case -> bailout
		marker <- 5
		fits <- c(rep(1,4), -check, rep(0, 4 + length(object@criteria)))
		return(list(fits = fits, marker = marker, done = done))		
	}
	diag(Phi) <-1
	object@Phi <- Phi

	check <- eigen(Phi, TRUE, TRUE)$values[ncol(Phi)]
	if(check < lower) { # Phi not posdef -> bailout
		marker <- 6
		fits <- c(rep(1,5), check, rep(0, 3 + length(object@criteria)))
		return(list(fits = fits, marker = marker, done = done))		

	}

	# Make beta
	beta <- object@beta$beta
	beta[object@beta$free] <- par[object@beta$select]

	if(SEFA) { # enforce and check Howe (1955) conditions
		object@beta$fix_beta_args$coefs  <- beta
		object@beta$fix_beta_args$cormat <- object@Phi
		beta <- do.call(FAiR_fix_coefficients, args = object@beta$fix_beta_args)
		check <- FAiR_check_coefficients(beta, threshold = lower)
		if(check != 1) { # fail rank check -> bailout
			marker <- 7
			fits <- c(rep(1,6), check, rep(0, 2 + length(object@criteria)))
			return(list(fits = fits, marker = marker, done = done))
		}
	}
	object@beta$beta <- beta

	signs <- ifelse(colSums(beta) >= 0, 1, -1)
	if(any(signs != 1) && !done) { # factors not all positive -> bailout
		marker <- 8
		fits  <- c(rep(1,7), mean(signs), rep(0, 1 + length(object@criteria)))
		return(list(fits = fits, marker = marker, done = done))
	}
	else { # get lexical criteria
		fits <- c(rep(1,8), FAiR_lexical_driver(par, object, S, 
							Delta, Xi, lower = lower))
		marker <- which(fits != 1)[1]
	}
	
	if( (marker == length(fits)) && (object@model == "SEFA") ) {
		# Mark which elements of Delta and beta got squashed to zero
		squashed_Delta <- object@Delta$select
		squashed_Delta[squashed_Delta] <- c(!Delta) & object@Delta$free

		squashed_beta <- object@beta$select
		squashed_beta[squashed_beta] <- c(!beta) & object@beta$free

		squashed <- squashed_Delta | squashed_beta
		return(list(fits = fits, marker = marker, 
			squashed = squashed, done = done))
	}
	else return(list(fits = fits, marker = marker, done = done))
})

setMethod("bfgs_fitS4", signature(par = "numeric", object = "restrictions.2ndorder", 
	helper = "list"), definition = 
function(par, object, helper, S, lower) {
	SEFA <- object@model == "SEFA"
	if(helper$marker == length(helper$fits)) {
		if(SEFA) { # squash based on helper and pretend it's a CFA
			object@model <- "CFA"
			small <- par[helper$squashed]
			par[helper$squashed] <- 0
		}
		fits <- fitS4(par, object, S, lower = lower) # get lexical criteria
		if(helper$done) return(fits[length(fits)])
		marker <- which(fits != 1)[1]
		if(marker < length(helper$fits)) return(-.Machine$integer.max)

		if(SEFA) out <- fits[marker] - crossprod(small) # penalize if small != 0
		else     out <- fits[marker]
		return(out)
	}
	else {
		fits <- fitS4(par, object, S, lower = lower) # get lexical criteria
		marker <- which(fits != 1)[1]
		if(marker > helper$marker) return(.Machine$integer.max)
		return(fits[marker])
	}
})

setMethod("gr_fitS4", signature(par = "numeric", object = "restrictions.2ndorder",
	helper = "list"), definition = 
function(par, object, helper, S, lower) {
	SEFA <- object@model == "SEFA"
	if( (helper$marker == length(helper$fits)) && (object@method == "MLE") ) {
		if(SEFA) { # squash based on helper and pretend it's a CFA
			small <- par[helper$squashed]
			par[helper$squashed] <- 0
		}

		# Make Theta2
		uniquenesses <- par[object@Theta2$select]
		uniquenesses_min <- min(uniquenesses)
		if(uniquenesses_min < 0) { # Heywood case -> bailout
			return(rep(.Machine$double.eps, length(par)))
		}
		else diag(object@Theta2$Theta2) <- uniquenesses

		# Make Xi
		mark_start <- 1
		mark_end   <- 0.5 * object@factors[2] * (object@factors[2] - 1)
		Xi <- object@Xi
		Xi[upper.tri(Xi)] <- par[mark_start:mark_end]
		Xi <- Xi + t(Xi)
		Xi <- FAiR_nearPD(Xi, posd.tol = lower)
		object@Xi <- Xi
		if(attributes(Xi)$ev < lower) {
			return(rep(.Machine$double.eps, length(par)))
		}
	
		# Make Delta
		Delta <- object@Delta$Delta
		Delta[object@Delta$free] <- par[object@Delta$select]
	
		signs <- ifelse(colSums(Delta) >= 0, 1, -1)
		if(any(signs != 1) && !helper$done) { # 2nd order factors not all positive -> reflect
			return(rep(.Machine$double.eps, length(par)))
		}
	
		# Make Phi
		Phi <- crossprod(chol(Xi) %*% t(Delta))
		check <- max(diag(Phi))
		if(check > 1) { # Heywood case -> bailout
			return(rep(.Machine$double.eps, length(par)))
		}
		diag(Phi) <-1
	
		check <- eigen(Phi, TRUE, TRUE)$values[ncol(Phi)]
		if(check < lower) { # Phi not posdef -> bailout
			return(rep(.Machine$double.eps, length(par)))
		}
		object@Phi <- Phi

		# Make beta
		beta <- object@beta$beta
		beta[object@beta$free] <- par[object@beta$select]

		signs <- ifelse(colSums(beta) >= 0, 1, -1)
		if(any(signs != 1) && !helper$done) { # factors not positive -> bailout
			return(rep(.Machine$double.eps, length(par)))
		}
		object@beta$beta <- beta

		fits <- FAiR_lexical_driver(par, object, S, Delta, Xi, lower = lower)
		marker <- which(fits != 1)[1]
		if(marker < length(fits)) { # failed a constraint -> bailout 
			return(rep(.Machine$double.eps, length(par)))
		}

		gradient <- FAiR_gradient(cormat = S, do_Phi = TRUE, 
						Phi = Phi, beta = beta, 
						Theta2 = object@Theta2$Theta2)
		dF_dPhi <- gradient$dF_dPhi
		diag(dF_dPhi) <- 0
		matderiv <- t(2 * Xi %*% t(Delta))
		dF_dDelta <- matrix(NA_real_, nrow = nrow(Delta), ncol = ncol(Delta))
		dF_dXi <- matrix(NA_real_, nrow(Xi), ncol(Xi))
		for(i in 1:ncol(Phi)) {
			dF_dDelta[i,] <- colSums(sweep(matderiv, 1, dF_dPhi[,i], 
							FUN = "*", check.margin = FALSE))
		}
		uppers <- which(upper.tri(Xi), arr.ind = TRUE)
		for(k in 1:nrow(uppers)) {
			d <- 0
			row <- uppers[k,1]
			col <- uppers[k,2]
			for(i in 1:nrow(Delta)) for(j in 1:nrow(Delta)) {
				d <- d + dF_dPhi[i,j] * Delta[i,row] * Delta[j,col] +
					 dF_dPhi[j,i] * Delta[i,col] * Delta[j,row]
			}
			dF_dXi[row,col] <- dF_dXi[col,row] <- d
		}

		gradient <- c(dF_dXi[upper.tri(Xi)], dF_dDelta[object@Delta$free], 
				gradient$dF_dbeta[object@beta$free], gradient$dF_dTheta2)
		if(SEFA) { # quadratic loss around zero for small coefficients if !done
			gradient[helper$squashed] <- -2 * small * (!helper$done)
		}
	}
	else { # find numeric gradient of whatever the marginal criterion is
		gradient <- FAiR_numeric_gradient(par, object, helper, S)
	}
	return(gradient)
})

setMethod("create_start", signature(restrictions = "restrictions.2ndorder"),
definition = function(restrictions, number, start, S) {
	factors  <- restrictions@factors[1]
	factors2 <- restrictions@factors[2]
	if(restrictions@model == "CFA") {
		out <- apply(restrictions@Domains, 1, FUN = function(x) {
				return(runif(number, min = max(x[1], -0.75), 
						     max = min(x[2],  0.75)))
			})
		return(out)		
	}
	Lambda  <- FAiR_get_loadings(1 - c(start), S, factors)
	signs <- sign(colSums(Lambda))
	if(any(signs != 1)) Lambda <- sweep(Lambda, 2, signs, FUN = "*", FALSE)
	Lambda <- Lambda %*% t(solve(FAiR_Landahl(factors)))
	out  <- matrix(NA_real_, nrow = number, ncol = restrictions@nvars)
	mark1 <- 0.5 * factors2 * (factors2 - 1)
	mark2 <- mark1 + sum(restrictions@Delta$free)
	while(number) {
		draws <- runif(mark2, min = restrictions@Domains[1:mark2, 1],
				      max = restrictions@Domains[1:mark2, 2])
		draws <- draws / 2
		Xi <- diag(factors2)
		Xi[upper.tri(Xi)] <- draws[1:mark1]
		Xi <- Xi + t(Xi)
		diag(Xi) <- 1
		Xi <- FAiR_nearPD(Xi)
		Delta <- restrictions@Delta$Delta
		Delta[restrictions@Delta$free] <- draws[(mark1 + 1):mark2]
		signs <- sign(colSums(Delta))
		if(any(signs != 1)) Delta <- sweep(Delta, 2, signs, FUN = "*", TRUE)
		FC <- Delta * (Delta %*% Xi)
		sorter2 <- order(colSums(FC), decreasing = TRUE)
		Phi <- try(crossprod(chol(Xi) %*% t(Delta)), TRUE)
		if(!is.matrix(Phi)) next
		if(any(diag(Phi) > 1)) next
		diag(Phi) <- 1
		Tmat <- try(chol(Phi), silent = TRUE)
		if(!is.matrix(Tmat)) next
		beta <- Lambda %*% t(chol2inv(Tmat))
		signs <- sign(colSums(beta))
		if(any(signs != 1)) beta <- sweep(beta, 2, signs, FUN = "*", FALSE)
		FC <- beta * (beta %*% Phi)
		sorter <- order(colSums(FC), decreasing = TRUE)
		beta <- beta[,sorter]
		Delta <- Delta[sorter,sorter2]
		Xi <- Xi[sorter2, sorter2]
		uniquenesses <- rnorm(nrow(start), 1 - start, sd = .1)
		out[number,] <- c(Xi[upper.tri(Xi)], Delta[restrictions@Delta$free], 
				  beta[restrictions@beta$free], uniquenesses)
		out[number,] <- ifelse(out[number,] < restrictions@Domains[,1], 
						      restrictions@Domains[,1],
				ifelse(out[number,] > restrictions@Domains[,2], 
						      restrictions@Domains[,2], 
						      out[number,]))
		number <- number - 1
	}
	return(out)
})

setMethod("create_FAobject", signature(restrictions = "restrictions.2ndorder"), 
definition = function(restrictions, opt, manifest, call, scores, lower) {
	S <- manifest$cor
	par <- opt$par
	factors <- restrictions@factors[1]
	factors2 <- restrictions@factors[2]
	SEFA <- restrictions@model == "SEFA"
	fits <- opt$value
	marker <- which(fits != 1)[1]
	if(marker < length(fits)) {
		warning("constraints did not bind")
	}

	# Make Theta2
	uniquenesses <- par[restrictions@Theta2$select]

	# Make Xi
	mark_start <- 1
	mark_end   <- 0.5 * restrictions@factors[2] * (restrictions@factors[2] - 1)
	Xi <- restrictions@Xi
	Xi[upper.tri(Xi)] <- par[mark_start:mark_end]
	Xi <- Xi + t(Xi)
	Xi <- FAiR_nearPD(Xi, posd.tol = lower)

	# Make Delta
	Delta <- restrictions@Delta$Delta
	Delta[restrictions@Delta$free] <- par[restrictions@Delta$select]
	if(SEFA) { # enforce and check Howe (1955) conditions on Delta
		restrictions@Delta$fix_Delta_args$coefs <- Delta
		restrictions@Delta$fix_Delta_args$cormat <- Xi
		Delta <- do.call(FAiR_fix_coefficients, args =restrictions@Delta$fix_Delta_args)
	}

	# Make Phi
	Phi <- crossprod(chol(Xi) %*% t(Delta))
	diag(Phi) <-1

	# Make beta
	beta <- restrictions@beta$beta
	beta[restrictions@beta$free] <- par[restrictions@beta$select]

	if(SEFA) { # enforce and check Howe (1955) conditions on beta
		restrictions@beta$fix_beta_args$coefs  <- beta
		restrictions@beta$fix_beta_args$cormat <- Phi
		beta <- do.call(FAiR_fix_coefficients, args = restrictions@beta$fix_beta_args)
	}

	# Make loadings_2nd and correlations_2nd
	Pi2 <- Delta %*% Xi
	FC2 <- Pi2 * Delta
	sorter_2nd <- order(colSums(FC2), decreasing = TRUE)
	uniquenesses_2nd <- 1 - rowSums(FC2)

	Xi_inv <- chol2inv(chol(Xi))
	D2 <- 1/sqrt(diag(Xi_inv))
	Psi2 <- sweep(sweep(Xi_inv, 1, D2, FUN = "*"), 2, D2, FUN = "*")
	Upsilon2 <- sweep(Delta, 2, D2, FUN = "*")
	RP2 <- sweep(Pi2, 2, D2, FUN = "/")

	loadings_2nd <- array(cbind(Delta, Pi2, RP2, Upsilon2, FC2),
				dim = c(factors, factors2, 5),
				dimnames = list(NULL, NULL, 
					c("PP", "PS", "RP", "RS", "FC")))

	correlations_2nd <- array(cbind(Xi, Psi2, diag(D2)),
			dim = c(factors2, factors2, 3), 
			dimnames = list(NULL, NULL, c("PF", "RF", "PR")))

	# construct primary pattern and factor contributions
	Pi <- beta %*% Phi
	FC <- beta * Pi
	uniquenesses <- 1 - rowSums(FC)
	names(uniquenesses) <- rownames(S)
	sorter <- order(colSums(FC), decreasing = TRUE)

	# construct reference structure and reference pattern
	Phi_inv <- chol2inv(chol(Phi))
	D <- 1/sqrt(diag(Phi_inv))
	Psi <- sweep(sweep(Phi_inv, 1, D, FUN = "*"), 2, D, FUN = "*")
	Upsilon <- sweep(beta, 2, D, FUN = "*")
	RP <- sweep(Pi, 2, D, FUN = "/")

	loadings <- array(cbind(beta, Pi, RP, Upsilon, FC),
			dim = c(nrow(S), factors, 5),
			dimnames = list(rownames(S), NULL, 
					c("PP", "PS", "RP", "RS", "FC")))

	correlations <- array(cbind(Phi, Psi, diag(D)),
			dim = c(factors, factors, 3), 
			dimnames = list(NULL, NULL, c("PF", "RF", "PR")))

	trans_mats <- array(diag(factors), dim = c(factors, factors, 2), 
			dimnames = list(paste("F", 1:factors), paste("F", 1:factors, 
					sep = ""), c("primary", "reference")))

	if(is.numeric(manifest$n.obs)) {
		par <- sweep(loadings[,,1], 1, sqrt(diag(manifest$cov)), FUN = "*", FALSE)
		par <- c(Xi[upper.tri(Xi)], Delta[restrictions@Delta$free], 
			par[restrictions@beta$free], uniquenesses * diag(manifest$cov))
		free <- par != 0
		Hessian <- FAiR_Hessian(par, restrictions, manifest$cov, lower)
		negHessian <- -1 * Hessian[free,free]
		ev <- eigen(negHessian, TRUE, TRUE)$values
		if(ev[length(ev)] < 0) {
			warning("Hessian indefinite, trying BHHH estimator")
			vcov   <- FAiR_BHHH(par, restrictions, manifest)
		}
		else {
			vcov <- try(chol2inv(chol(negHessian)))
			if(is.matrix(vcov)) vcov <- vcov * (2/(manifest$n.obs))
			else vcov <- FAiR_BHHH(par, restrictions, manifest)
		}
		zstats  <- par[free] / sqrt(diag(vcov))
		zstats_Xi <- NA_real_ * Xi
		mark_start <- 1
		mark_end <- 0.5 * factors2 * (factors2 - 1)
		zstats_Xi[upper.tri(zstats_Xi)] <- zstats[mark_start:mark_end]
		zstats_Xi <- zstats_Xi + t(zstats_Xi)
		diag(zstats_Xi) <- NA_real_
		zstats_Xi <- zstats_Xi[sorter_2nd, sorter_2nd]
		fill <- restrictions@Delta$free & c(Delta != 0)
		mark_start <- mark_end + 1
		mark_end <- mark_end + sum(fill)
		zstats_Delta <- NA_real_ * Delta
		zstats_Delta[fill] <- zstats[mark_start:mark_end]
		zstats_Delta <- zstats_Delta[sorter,sorter_2nd]
		fill <- restrictions@beta$free & c(beta != 0)
		zstats_beta <- NA_real_ * beta
		mark_start <- mark_end + 1
		mark_end <- mark_end + sum(fill)
		zstats_beta[fill] <- zstats[mark_start:mark_end]
		zstats_beta <- zstats_beta[,sorter]
		mark_start <- mark_end + 1
		mark_end <- length(zstats)
		zstats_Theta2 <- zstats[mark_start:mark_end]
		zstats <- list(Xi = zstats_Xi, Delta = zstats_Delta, 
				beta = zstats_beta, Theta2 = zstats_Theta2)
	}
	else {
		warning("zstatistics could not be calculated because the number of",
			" observations is unknown")
		vcov <- matrix(NA_real_,nrow = restrictions@nvars, 
					ncol = restrictions@nvars)
		zstats <- list()
	}

	loadings_2nd <- loadings_2nd[sorter,sorter_2nd,]
	correlations_2nd <- correlations_2nd[sorter_2nd,sorter_2nd,]
	loadings <- loadings[,sorter,]
	correlations <- correlations[sorter,sorter,]
	attributes(correlations)$orthogonal <- FALSE
	rownames(correlations) <- colnames(correlations) <- colnames(loadings) <- 
							paste("F", 1:factors, sep = "")

	scores <- FAiR_scores(scores, manifest$zz, beta, Phi, uniquenesses, S)
	
	if(is.null(call$seeds)) seeds <- rep(formals(Factanal)$seeds, 2)
	else {
		seeds <- eval(call$seeds)
		if(length(seeds) == 1) seeds <- c(seeds, seeds)
	}

	seeds <- matrix(seeds, nrow = 1)
	colnames(seeds) <- c("unif.seed", "int.seed")

	loadings_2nd <- loadings_2nd[sorter,,]
	if(SEFA) { # recalculate degrees of freedom
		dof <- 0.5 * ncol(S) * (ncol(S) + 1) - restrictions@nvars
		dof <- dof + sum(restrictions@Domains[,1] == restrictions@Domains[,2])
		dof <- dof + sum(!beta) + sum(!Delta)
		dof <- dof + sum(restrictions@beta$beta   != 0,   na.rm = TRUE)
		dof <- dof + sum(restrictions@Delta$Delta != 0, na.rm = TRUE)
		restrictions@dof <- as.integer(dof)
	}

	FAobject <- new("FA.2ndorder", loadings = loadings, 
			loadings_2nd = loadings_2nd, correlations = correlations,
			correlations_2nd = correlations_2nd, 
			trans_mats = trans_mats, uniquenesses = uniquenesses,
			uniquenesses_2nd = uniquenesses_2nd,
			restrictions = restrictions, vcov = vcov, zstats = zstats,
			scores = scores, manifest = manifest, rgenoud = list(extraction = opt), 
			model = restrictions@model, method = restrictions@method, 
			call = call,  seeds = seeds)
	return(FAobject)
})


## The rest are various common methods; logLik and summary are established in stats4
setMethod("coef", "FA", function(object, matrix = "PP", ...) { 
	as.matrix(object@loadings[,,toupper(matrix)])
})
setMethod("coef", "FA.general", function(object, matrix = "PP", level = 1, ...) {
	if(level == 1) return(as.matrix(object@loadings[,,toupper(matrix)]))
	else           return(as.matrix(object@loadings_2nd))
})

setMethod("coef", "FA.2ndorder", function(object, matrix = "PP", level = 1, ...) {
	if(level == 1) return(as.matrix(object@loadings[,,toupper(matrix)]))
	else           return(as.matrix(object@loadings_2nd[,,toupper(matrix)]))
})

setMethod("vcov", "FA", function(object, ...) object@vcov)

setMethod("logLik", "FA",
function (object, ...) {
	if(object@method != "MLE") stop("log-likelihood can only be calculated for MLEs")
        LHS <- model.matrix(object)
	RHS <- fitted(object) + diag(object@uniquenesses)
	RHS_inverse <- chol2inv(chol(RHS))

	p <- ncol(LHS)
        n.obs <- object@manifest$n.obs

	val <-  -0.5 * (n.obs - 1) * ( log(det(RHS)) + sum(diag(LHS %*% RHS_inverse)) )
	
	if(object@model == "EFA") {
		Lambda <- coef(object)
		npars <- length(Lambda) + nrow(Lambda)
		correction <- 0.5 * ncol(Lambda) * (ncol(Lambda) - 1)
	}
	else if(object@model == "SEFA") {
		npars <- object@restrictions@nvars
		correction <- sum( c(!coef(object)) & object@restrictions@beta$free )
		if(is(object, "FA.2ndorder")) {
			correction <- correction + sum( c(!coef(object, level = 2)) &
						object@restrictions@Delta$free )
		}
	}
	else { # CFA
		npars <- object@restrictions@nvars
		correction <- 0
	}
	attr(val, "df")   <- npars - correction
	attr(val, "nobs") <- n.obs
	class(val)        <- "logLik"
	return(val)
})

setMethod("BIC", signature(object = "FA"),
function (object, ...) {
        BIC(logLik(object))
})

setMethod("confint", signature(object = "FA"), definition = 
function (object, parm, level = 0.95, ...){
        stop("confint does not work yet")
})

setMethod("profile", signature(fitted = "FA"), definition = 
function (fitted, delta = 0.05, number = 100, plot.it = TRUE, ...) {
	if(attributes(fitted@correlations)$orthogonal) {
		stop("profile is not defined for orthogonal solutions")
	}
	if(fitted@method != "MLE") {
		stop("profile is currently only defined for maximum-likeliood models")
	}
	notfree <- !fitted@restrictions@beta$free
	notfree[!coef(fitted)] <- TRUE
	notfree <- matrix(notfree, ncol = fitted@restrictions@factors[1])
	width <- seq(from = -delta, to = delta, length = number)
	outlist <- list()
	count <- 1
	Theta2 <- diag(fitted@uniquenesses)
	parnames <- rownames(coef(fitted))
	ask <- par()$ask
	par("ask" = TRUE)
	on.exit(par("ask" = ask))
	for(p in 1:ncol(notfree)) for(j in 1:nrow(notfree)) if(notfree[j,p]) {
		y <- x <- rep(NA_real_, length(width))
		val <- fitted@loadings[j,p,1]
		for(i in 1:length(x)) {
			x[i] <- fitted@loadings[j,p,1] <- val + width[i]
			y[i] <- logLik(fitted)[1]
		}
		outlist[[count]] <- list(x = x, y = y)
		names(outlist)[count] <- paste(parnames[j], p, sep = "_")
		fitted@loadings[j,p,1] <- val
		if(plot.it) {
			plot(x, y, type = "l", ylab = "Log-likelihood",
				xlab = paste("Factor", p, "for", parnames[j]), ...)
			abline(v = val, col = "gray", lty = "dotted")
		}
		count <- count + 1
	}
	return(invisible(outlist))
})

setMethod("show", "FA", 
function(object) {
	cat("\nCall:\n")
	print(object@call)
	show(object@restrictions)
	x <- FAiR_stats(object)
	cat("\nLog-likelihood for model: ", x$llik, "\n")
	cat("p-value for null hypothesis that model holds in the population: ", x$p.value)
        cat("\nSIC for model:", x$SIC)
	cat("\nBIC for model: ", x$BIC)
	cat("\nBIC for a saturated model: ", x$BIC_saturated)
	## The following is slightly modified from sem:::print.summary.sem, which is
	## Copyright 2007 John Fox and is licensed under the GPL V2+
	cat("\nGoodness-of-fit index = ", x$GFI)
	cat("\nAdjusted goodness-of-fit index = ", x$AGFI)
        cat("\nRMSEA index =  ", x$RMSEA[1],
            "   ", 100*x$RMSEA[4], "% CI: (", x$RMSEA[2], ", ", x$RMSEA[3],")", sep="")
        cat("\nBentler-Bonnett NFI = ", x$NFI)
        cat("\nTucker-Lewis NNFI = ", x$NNFI)
        cat("\nBentler CFI = ", x$CFI)
	cat("\n")
})

setMethod("plot", "FA", 
function(x, y, ...) {
	if(!require(nFactors)) {
		stop("plot requires that the nFactors package to be installed")
	}
	S <- model.matrix(x)
	eigenvalues <- eigen(S, TRUE, TRUE)$values
	variables <- length(eigenvalues)
	nsubjects <- x@manifest$n.obs
	if(!is.numeric(nsubjects)) {
		stop("parallel analysis is only possible when the number of observations is known\n",
			"please call Factanal() again passing the raw data through 'x' or\n",
			"passing the number of observations through 'n.obs'", call. = FALSE)
	}
	aparallel <- parallel(var = variables, subject = nsubjects)$eigen$qevpea
	results   <- nScree(eig = eigenvalues, aparallel = aparallel)
	plotnScree(results)
	return(invisible(results))
})

setMethod("summary", signature("FA"), def =
function(object, ...){
	if(all(diag(object@manifest$cov) == 1)) {
		warning("z-statistics cannot be accurately calculated yet when a correlation matrix ",
			" is modeled\nbecause the standard errors do not take into account the fact",
			" that the sample\nstandard deviations are estimates. For better",
			" z-statistics, reestimate the model after \npassing either the raw data or",
			" the sample covariance matrix to Factanal()")
	}
	if(attributes(object@correlations)$orthogonal) {
		return(new("summary.FA", point_estimates = list(object@uniquenesses),
					 zstats = list(object@zstats[["Theta2"]]),
					 restrictions = object@restrictions,
					 call = object@call))
	}
	else {
		if(is(object@restrictions, "restrictions.factanal")) {
			point_estimates <- list(beta = coef(object),
						Phi  = as.matrix(object@correlations[,,"PF"]),
						Theta2 = object@uniquenesses)
			zstats <- list(beta = object@zstats$beta,
					Phi = object@zstats$Phi,
					Theta2 = object@zstats$Theta2)
		}
		else if(is(object@restrictions, "restrictions.orthonormal")) {
			point_estimates <- list(beta = coef(object),
						Phi  = object@correlations[,,"PF"],
						Theta2 = object@uniquenesses)
			zstats <- list(beta = object@zstats$beta,
					Phi = object@zstats$Phi,
					Theta2 = object@zstats$Theta2)
		}
		else if(is(object@restrictions, "restrictions.1storder")) {
			point_estimates <- list(beta = coef(object),
						Phi  = object@correlations[,,"PF"],
						Theta2 = object@uniquenesses)
			zstats <- list(beta = object@zstats$beta,
					Phi = object@zstats$Phi,
					Theta2 = object@zstats$Theta2)
		}
		else if(is(object@restrictions, "restrictions.general")) {
			point_estimates <- list(beta = coef(object),
						Phi  = object@correlations[,,"PF"],
						Theta2 = object@uniquenesses,
						Delta = coef(object, level = 2))
			zstats <- list(beta = object@zstats$beta,
					Delta = object@zstats$Delta,
					Theta2 = object@zstats$Theta2)
		}
		else if(is(object@restrictions, "restrictions.2ndorder")) {
			point_estimates <- list(beta = coef(object),
						Phi  = object@correlations[,,"PF"],
						Theta2 = object@uniquenesses,
						Delta = coef(object, level = 2),
						Xi = object@correlations_2nd[,,"PF"])
			zstats <- list(beta = object@zstats$beta,
					Delta = object@zstats$Delta,
					Theta2 = object@zstats$Theta2,
					Xi = object@zstats$Xi)
		}
		else stop("no summary method available")

		return(new("summary.FA", point_estimates = point_estimates,
					 zstats          = zstats,
					 restrictions = object@restrictions,
					 call = object@call))
	}
})

setMethod("show", "summary.FA", function(object){
	if(length(object@point_estimates) == 1) {
		cat("\nCall:\n")
		print(object@call)
		mat <- as.matrix(object@point_estimates[[1]])
		if(!is.null(object@zstats[[1]])) { 
			mat <- cbind(mat, object@zstats[[1]])
			colnames(mat) <- c("uniquenesses", "z-statistics")
		}
		else colnames(mat) <- "uniquenesses"
		print(round(mat, 3))
	}
	else {
		point_estimates <- cbind(object@point_estimates$beta, 0, 
					 object@point_estimates$Theta2)
		if(ncol(object@point_estimates$Phi) > 1) {
			point_estimates <- rbind(point_estimates, 0,
				   		cbind(object@point_estimates$Phi, 0, 0))
		}

		if(!is.null(object@point_estimates$Xi)) {
			point_estimates <- rbind(point_estimates, 0, cbind(
						t(object@point_estimates$Delta), 0, 0))
			point_estimates <- rbind(point_estimates, 0, cbind(
						object@point_estimates$Xi, matrix(0,  
						nrow = nrow(object@point_estimates$Xi),
						ncol = ncol(point_estimates) - 
						       ncol(object@point_estimates$Xi))))
			if(!all(sapply(object@zstats, FUN = is.null))) {
				zstats <- cbind(object@zstats$beta, 0, object@zstats$Theta2)
				zstats <- rbind(zstats, 0, cbind(t(object@zstats$Delta), 0, 0))
				zstats <- rbind(zstats, 0, cbind(object@zstats$Xi, matrix(0,
						nrow = nrow(object@zstats$Xi), ncol = 
						ncol(zstats) - ncol(object@zstats$Xi))))
				rownames(zstats)[which(rownames(point_estimates) == "2nd order")] <-
										    "2nd order"
			}
			else zstats <- NULL
		}
		else if(!is.null(object@point_estimates$Delta)) {
			point_estimates <- rbind(point_estimates, 0, cbind(
						t(object@point_estimates$Delta), 0, 0))
			rownames(point_estimates)[nrow(point_estimates)] <- "2nd order"
			if(!all(sapply(object@zstats, FUN = is.null))) {
				zstats <- cbind(object@zstats$beta, 0, object@zstats$Theta2)
				zstats <- rbind(zstats, 0, cbind(t(object@zstats$Delta), 0, 0))
				rownames(zstats)[nrow(zstats)] <- "2nd order"
			}
			else zstats <- NULL
		}
		else if(!is.null(object@zstats$Phi)) {
			zstats <- cbind(object@zstats$beta, 0, object@zstats$Theta2)
			zstats <- rbind(zstats, 0, cbind(object@zstats$Phi, 0, 0))
			rownames(zstats) <- rownames(point_estimates)
		}
		else zstats <- NULL

		colnames(point_estimates) <- c(paste("F", 
						1:ncol(object@point_estimates$beta), 
						sep = ""), "", "Uniqueness")

		## This part is slightly modified from stats:::print.loadings, which is 
        	## Copyright 1995-2007 R Core Development Team and licensed under GPL V2+
		cat("\nCall:\n")
		print(object@call)
		cat("\nPoint estimates (blanks, if any, are exact zeros):\n")
		fx <- format(round(point_estimates, 3))
		names(fx) <- NULL
		nc <- nchar(fx[1], type = "c")
		fx[!point_estimates] <- paste(rep(" ", nc), collapse = "")
		print(fx, quote = FALSE)
		if(!is.null(zstats)) {
			colnames(zstats) <- colnames(point_estimates)
			cat("\nz-statistics for null hypothesis that estimate is zero:\n")
			fx <- format(round(zstats, 3))
			names(fx) <- NULL
			nc <- nchar(fx[1], type = "c")
			fx[!zstats] <- paste(rep(" ", nc), collapse = "")
			print(fx, quote = FALSE)
		}
	}
})

setMethod("plot", signature(x = "summary.FA"), def = 
function(x, y = NULL, ...) {
	## Do something useful, maybe ggplot() or whatever the new thing is
})

setMethod("show", signature(object = "restrictions"), definition = 
function (object) {
	str(object)
})

setMethod("show", signature(object = "restrictions.factanal"), definition = 
function (object) {
	cat("\nExploratory factor anaylsis\n")
	cat(object@factors[1], "factors\n")
	cat(object@dof, "degrees of freedom\n")
})

setMethod("show", signature(object = "restrictions.orthonormal"), definition = 
function (object) {
	cat("\nExploratory factor anaylsis\n")
	cat(object@factors[1], "factors\n")
	cat(object@dof, "degrees of freedom\n")
})

setMethod("show", signature(object = "restrictions.1storder"), definition = 
function (object) {
	if(object@model == "SEFA") {
		cat("\nSemi-exploratory factor anaylsis\n")
		FAiR_show_coef(object, level = 1)
		FAiR_show_sefa(object@beta$fix_beta_args, level = 1)
	}
	else {
		cat("\nConfirmatory factor anaylsis\n")
		FAiR_show_coef(object, level = 1)	
	}
	factors <- object@factors[1]
	mark_start <- 1
	mark_end <- factors
	Phi_Domains <- object@Domains[mark_start:mark_end,,drop = FALSE]
	if(any(Phi_Domains[,1] > (-1 + .Machine$double.eps)) | 
	   any(Phi_Domains[,1] > ( 1 - .Machine$double.eps)) ) {
		cat("\nBounds on correlations among primary factors\n")
		print(Phi_Domains)
	}
	else cat("\nAll correlations among primary factors in the [-1,1] interval\n")
	FAiR_show_constraints(object@criteria)
	cat("\n", object@dof, "degrees of freedom\n")
})

setMethod("show", signature(object = "restrictions.general"), definition = 
function (object) {
	if(object@model == "SEFA") {
		cat("\nSemi-exploratory factor anaylsis\n")
		FAiR_show_coef(object, level = 1)
		FAiR_show_sefa(object@beta$fix_beta_args, level = 1)
		FAiR_show_coef(object, level = 2)
	}
	else { # CFA
		cat("\nConfirmatory factor anaylsis\n")
		FAiR_show_coef(object, level = 1)
		FAiR_show_coef(object, level = 2)
	}
	FAiR_show_constraints(object@criteria)
	cat("\n", object@dof, "degrees of freedom\n")
})

setMethod("show", signature(object = "restrictions.2ndorder"), definition = 
function (object) {
	if(object@model == "SEFA") {
		cat("\nSemi-exploratory factor anaylsis\n")
		FAiR_show_coef(object, level = 1)
		FAiR_show_sefa(object@beta$fix_beta_args, level = 1)
		FAiR_show_coef(object, level = 2)
	}
	else { # CFA
		cat("\nConfirmatory factor anaylsis\n")
		FAiR_show_coef(object, level = 1)
		FAiR_show_coef(object, level = 2)
	}
	factors <- object@factors
	mark_start <- 1
	mark_end <- factors[2]
	Xi_Domains <- object@Domains[mark_start:mark_end,,drop=FALSE]
	if(any(Xi_Domains[,1] > (-1 + .Machine$double.eps)) | 
	   any(Xi_Domains[,1] > ( 1 - .Machine$double.eps)) ) {
		cat("\nBounds on correlations among primary factors at level 2\n")
		print(Xi_Domains)
	}
	else cat("\nAll correlations among primary factors at level 2 ",
		"in the [-1,1] interval\n")
	FAiR_show_constraints(object@criteria)
	cat("\n", object@dof, "degrees of freedom\n")
})
