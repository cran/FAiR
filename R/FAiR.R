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

## This file defines functions that are meant to be accessed by users directly.

## NOTE: This file is intended to be read with 90 columns and 8 space tabs

## Factanal() estimates all factor analysis models
Factanal <-
function(x, factors, data = NULL, covmat = NULL, n.obs = NA, subset, na.action, 
	scores = "none", seeds = 12345, lower = sqrt(.Machine$double.eps), 
	model = c("SEFA", "EFA", "CFA"), method = c("MLE", "YWLS"), 
	restrictions, fixed, criteria = NULL, robust.covmat = FALSE, ...) {

## Arguments x through scores are basically the same as in stats:::factanal
#  x: a formula or matrix with outcome variables, can be ordinal in the case of a formula
#  factors: the number of factors (though it can now be a vector of length 2)
#  data: an optional dataframe if x is a formula
#  covmat: a covariance matrix or a list with an element named cov
#  n.obs: the number of observations if covmat is used to pass a matrix
#  subset: a logical vector indicating which observations to use among x 
#  na.action: character indicating what to do with if x has NA
#  scores: character indicating whether / how to calculate factor scores

#  seeds: PRNG seeds for the unif.seed and int.seed arguments of genoud()
#  lower: lower bound on uniquenesses in factanal()-style EFA models and the lower bound
#         on eigenvalues when determining whether a matrix is computationally posdef
#  method: character indicating what estimation method to use
#  model:  character indicating what model to estimate
#  restrictions: object of class "restrictions" (usually better to leave unspecified)
#  fixed: an optional matrix of values to fix certain coefficients in CFA or SEFA models
#  criteria: an optional character vector or list of characters indicating what criteria
#            to include under lexical optimization (usually better to leave unspecified)
#  robust: logical indicating whether the minimum covariance determinant method should
#          be used to estimate the covariance of the outcomes
#  robust.covmat = logical indicating whether to use minimum distance covariance estimator
#    ...: arguments that get passed to genoud()

	## Preliminaries
	stopifnot(require(rgenoud))

        kall <- match.call()
	manifest <- FAiR_parse(x = x, data = data, covmat = covmat, n.obs = n.obs,
				subset = subset, na.action = na.action,  
				robust = robust.covmat, seeds = seeds)
				
	S <- manifest$cor

	scores   <-  match.arg(scores, c("none", "regression", "Bartlett", "Thurstone",
					"Ledermann", "Anderson-Rubin", "McDonald",
					"Krinjen", "Takeuchi", "Harman"))

	if(missing(restrictions)) {
		if(factors[1] == 1) {
			if(model != "EFA") {
				warning("one factor was requested so the most efficent ",
					"EFA algorithm will be used to extract it")
				model  <- "EFA"
				method <- "MLE"
			}
			else if(method != "MLE") {
				warning("one factor was requested so MLE will be used ",
					"to extract it for computational efficiency")
				method <- "MLE"
			}
			p <- nrow(S)
			restrictions <- new("restrictions.factanal", 
						factors = c(1,0), nvars = p,
						dof = as.integer(0.5 * p * (p + 1) - 2 * p),
						Domains = cbind(lower, rep(1,p)),
						model = model, method = method)	
		}
		else {
			model    <-  match.arg(toupper(model),  eval(formals(Factanal)$model))
			method   <-  match.arg(toupper(method), eval(formals(Factanal)$method))
			restrictions <- make_restrictions(factors, model, method, 
							  fixed, S, criteria)
		}
	}
	else {
		if(!is(restrictions, "restrictions")) {
			stop("'restrictions' must be an object of class 'restrictions'",
				" or left unspecified (recommended)")
		}
		model <- restrictions@model
		method <- restrictions@method
		factors <- restrictions@factors
	}

	if(model == "EFA") {
		if(!missing(fixed)) {
			stop("providing 'fixed' is inconsistent with exploratory factor analysis")
		}
	}

	## Prepare to call genoud via the model.frame trick
        mc <- match.call(expand.dots = TRUE)
        mc[[1]] <- as.name("genoud")
        mc[names(formals(Factanal))] <- NULL

	# Arguments for genoud() that are logically required by Factanal()
	mc$nvars    <- restrictions@nvars
	mc$max      <- TRUE
	mc$hessian  <- FALSE # we do our own Hessian inside create_FAobject()
	mc$lexical  <- TRUE
	mc$Domains  <- restrictions@Domains
	mc$default.domains <- NULL
	mc$data.type.int  <- FALSE

	mc$fn <- function(par) {
			fitS4(par, restrictions, S, lower = lower)
		}
	mc$BFGSfn <- function(par, helper = NA) {
			bfgs_fitS4(par, restrictions, helper, S, lower = lower)
		}
	mc$BFGShelp <- function(initial, done = FALSE) {
				flush.console()
			 	bfgs_helpS4(initial, restrictions, done, S, lower = lower)
		}
	mc$gr <- if(method == "MLE") function(par, helper = NA) {
		 gr_fitS4(par, restrictions, helper, S, lower = lower)
		} else NULL
	
	if(is(restrictions, "restrictions.factanal")) { # Deal with some EFA stuff
		if(restrictions@fast) {
			efa <- factanal(x, factors[1], data, covmat, n.obs, subset, na.action,
					scores = "none", rotation = "none", 
					control = list(lower = lower))
			mc$starting.values <- matrix(efa$uniquenesses, nrow = 1)
			mc$gradient.check <- FALSE
			mc$max.generations <- 1
			mc$BFGSburnin <- 2
		}
		mc$BFGSfn   <- NULL
		mc$BFGShelp <- NULL
		mc$lexical  <- FALSE
		mc$boundary.enforcement <- 2
		mc$Domains[,1] <- pmax(restrictions@Domains[,1], lower)
	}

	# Workaround to get replicatability from genoud()
	if(any(!is.null(mc$unif.seed) | !is.null(mc$int.seed))) {
		warning("Use the seeds argument to Factanal() instead of the unif.seed",
			"and int.seed arguments to genoud(). Using 12345 as the seed.")
	}
	mc$unif.seed <- seeds[1]
	mc$int.seed  <- if(length(seeds) == 1) seeds else seeds[2]
	set.seed(mc$unif.seed)

	# Arguments for genoud() that are superceded unless explicitly specified
	if(is.null(mc$boundary.enforcement)) mc$boundary.enforcement <- 1
	if(is.null(mc$pop.size))             mc$pop.size <- formals(genoud)$pop.size
	if(is.null(mc$MemoryMatrix))         mc$MemoryMatrix <- FALSE
	if(is.null(mc$print.level))          mc$print.level <- 1
	if(is.null(mc$P9mix))                mc$P9mix <- 1
	if(is.null(mc$BFGSburnin))           mc$BFGSburnin <- -1
	if(is.null(mc$max.generations))      mc$max.generations <- 1000
	if(is.null(mc$project.path))         mc$project.path <- paste(tempfile(), 
								"Factanal.txt", sep = "")

	# Deal with starting values
	if(is.null(mc$starting.values)) {
		cat("Generating good starting values, have patience ...\n")
		flush.console()
		start <- as.matrix( FAiR_PACE_by_RGENOUD(S, factors[1], seeds = seeds) )
		pop.size <- eval(mc$pop.size)
		start <- create_start(restrictions, pop.size, start, S)

		if(pop.size >= nrow(start)) {
			mc$starting.values <- start[,1:eval(mc$nvars)]
		}
		else mc$starting.values <- start[nrow(start):(nrow(start) - pop.size + 1),
								1:mc$nvars]
					
	}
	else if(!is.matrix(eval(mc$starting.values))) {
		if(length(eval(mc$starting.values)) != ncol(S)) {
			stop("you need to supply a vector of length ", mc$nvars, " or ",  
			      ncol(S), " if starting.values is specified")
		}
		pop.size <- eval(mc$pop.size)
		start <- create_start(restrictions, pop.size, 
				      as.matrix(eval(mc$starting.values)), S)

		if(pop.size >= nrow(start)) {
			mc$starting.values <- start[,1:eval(mc$nvars)]
		}
		else mc$starting.values <- start[nrow(start):(nrow(start) - pop.size + 1),
								1:mc$nvars]		

	}
	else if(ncol(eval(mc$starting.values)) != mc$nvars) {
		if(ncol(eval(mc$starting.values)) != ncol(S)) {
			stop("you need to supply a matrix with ", mc$nvars, " or ",  
			      ncol(S), " columns if starting.values is specified")
		}
		pop.size <- eval(mc$pop.size)
		start <- create_start(restrictions, pop.size, eval(mc$starting.values), S)

		if(pop.size >= nrow(start)) {
			mc$starting.values <- start[,1:eval(mc$nvars)]
		}
		else mc$starting.values <- start[nrow(start):(nrow(start) - pop.size + 1),
								1:mc$nvars]		
	}
	opt <- eval(mc) # does all the real work

	## Call method to postprocess opt and bake object of correct class
	FAobject <- create_FAobject(restrictions, opt, manifest, kall, scores, lower)
	return(FAobject)
}

## Rotate() finds an oblique transformations of the factors following EFA extraction
Rotate <-
function(FAobject, criteria, weights = NULL, seeds = 12345, ...) {
## Arguments
# FAobject: the thing returned by Factanal()
# criteria: a list of functions or a list of character strings with names of functions or
#           leave it unspecified to get help from the GUI
# weights:  a vector of non-negative weights all less than or equal to 1.0 with length
#           equal to the number of factors less one to be used if the varphi criterion is
#           to be used as the lexical criterion in the lexical minimization. Leave as NULL
#           if varphi is not used or if varphi is used but dynamic weights are desired
#    seeds: PRNG seeds for the unif.seed and int.seed arguments of genoud()

	stopifnot(require(rgenoud))

	if(!is(FAobject, "FA")) {
		stop("FAobject must be produced by FA()")
	}
	if(FAobject@model != "EFA") {
 		stop("Rotate() only works for exploratory factor analysis")
	}

	if(missing(criteria)) {
		criteria <- make_criteria(ncol(coef(FAobject)), weights, FAobject)
	}
	else if(length(criteria) == 0) {
		criteria <- list(FAiR_criterion_varphi)
		formals(criteria[[1]])$weights <- weights
	}
	else for(i in 1:length(criteria)) {
		if(is.character(criteria[[i]])) {
			criteria[[i]] <- get(paste("FAiR_criterion_", 
						criteria[[i]], sep = ""))
		}
		else if(!is.function(criteria[[i]])) {
			stop("criteria must be a list of character strings or functions",
				"it is usually best to leave criteria unspecified")
		}
	}
	if(length(criteria) == 1) criteria <- c(list(FAiR_criterion_no_factor_collapse),
						criteria)

	## Prepare to call genoud via the model.frame trick
        mc <- match.call(expand.dots = TRUE)
        mc[[1]] <- as.name("genoud")
        mc[names(formals(Rotate))] <- NULL

	# Arguments for genoud() that are logically required by Rotate()
	mc$nvars    <- ncol(coef(FAobject))^2
	mc$max      <- FALSE
	mc$hessian  <- FALSE
	mc$lexical  <- TRUE
	mc$Domains  <- NULL
	mc$default.domains <- 1
	mc$data.type.int  <- FALSE

	mc$fn     <- function(par) FAiR_Rotate(par, coef(FAobject), criteria) 
	mc$BFGSfn <- function(par, helper) {
				fits <- FAiR_Rotate(par, coef(FAobject), criteria)
				marker <- which(fits != -1)[1]
				if(marker > helper$marker) return(-.Machine$double.xmax)
				if(marker < helper$marker) return( .Machine$double.xmax)
				                           return(fits[marker])
			}
 	mc$gr <- NULL
	mc$BFGShelp <- function(initial, done = FALSE) {
				fits <- FAiR_Rotate(initial, coef(FAobject), criteria)
				marker <- which(fits != -1)[1]
				return(list(fits = fits, marker = marker))
			}

	# Workaround to get replicatability from genoud()
	if(any(!is.null(mc$unif.seed) | !is.null(mc$int.seed))) {
		warning("Use the seeds argument to Factanal() instead of the unif.seed",
			"and int.seed arguments to genoud(). Using 12345 as the seed.")
	}
	mc$unif.seed <- seeds[1]
	mc$int.seed  <- if(length(seeds) == 1) seeds else seeds[2]
	set.seed(mc$unif.seed)

	# Arguments for genoud() that are superceded unless explicitly specified
        if(is.null(mc$pop.size))        mc$pop.size <- formals(genoud)$pop.size
	if(is.null(mc$print.level))     mc$print.level <- 1
	if(is.null(mc$max.generations)) mc$max.generations <- 1000
        if(is.null(mc$boundary.enforcement)) mc$boundary.enforcement <- 1
	if(is.null(mc$MemoryMatrix))    mc$MemoryMatrix <- FALSE
	if(is.null(mc$P9mix))           mc$P9mix <- 1
	if(is.null(mc$BFGSburnin))      mc$BFGSburnin <- -1
	if(is.null(mc$project.path))    mc$project.path <- paste(tempfile(), 
								"Rotate.txt", sep = "")
	if(is.null(mc$starting.values)) {
# 		starting.values <- FAiR_make_starts_Rotate(number = eval(mc$pop.size),
#                                                                 A = coef(FAobject))
		mc$starting.values <- NULL
	}
        else if(!is.matrix(eval(mc$starting.values))) {
		stop("starting.values must be a matrix or unspecified")
	}
	else if(ncol(eval(mc$starting.values)) != mc$nvars) {
		if(eval(nrow(mc$starting.values) == mc$nvars)) {
			mc$starting.values <- t(mc$starting.values)
		}
		else stop("starting.values has the wrong number of columns")
	}
	opt <- eval(mc)
	FAobject <- FAiR_opt2FAobject(opt, FAobject, seeds)
	return(FAobject)
}

## read.cefa() tries to read files produced by the Comprehensive Exploratory Factor
## Analysis (CEFA) package but is probably pretty fragile
read.cefa <- read.CEFA <- 
function(file) {
	lines <- readLines(file)

	n.obs <- as.integer(strsplit(lines[[1]], split = " ", extended = FALSE)[[1]][1])

	datatype <- grep("([D|d]ata", lines, extended = FALSE)
	if(length(datatype) == 0) datatype <- lines[2]
	datatype <- as.integer(strsplit(lines[datatype[1]], split = "", 
					extended = FALSE)[[1]][1])

	if(datatype == 3) {
		stop("importing a factor loading matrix is not supportable", 
			"please import the raw data or the covariance matrix")
	}

	varnames_marker <- grep("([V|v]ariable", lines, extended = FALSE)
	if(length(varnames_marker) == 0) varnames_marker <- 4
	else varnames_marker <- varnames_marker[1]

	varnames <- as.integer(strsplit(lines[varnames_marker], split = "",
					extended = FALSE)[[1]][1])
	if(varnames == 1) {
		varnames <- as.character(NULL)
		varnames_marker <- varnames_marker + 1
		while(varnames_marker) {
			if(length(grep("[a-z]", lines[varnames_marker], 
					ignore.case = TRUE, extended = FALSE)) > 0) {
					
					tempnames <- strsplit(lines[varnames_marker],
							split = " ", extended = FALSE)
					if(length(grep("^[\\s]*[0-9]", tempnames[[1]][1], 
						extended = FALSE)) > 0) break
					varnames <- c(varnames, unlist(sapply(tempnames,
							FUN = function(x) x[x!=""])))
					varnames_marker <- varnames_marker + 1
				}
			else if(length(grep("[0-9]", lines[varnames_marker],
					extended = FALSE)) > 0) break
			else {
				varnames <- NULL
				warning("there was problem assigning variable names")
				break
			}
		}
	}
	else if(varnames == 0) varnames <- NULL
	else {
		varnames <- NULL
		warning("could not determine whether there are variable names")
	}

	lines <- lines[-c(1:6)]
	mark <- grep("[a-z]", lines, ignore.case = TRUE, extended = TRUE)
	if(length(mark) > 0) lines <- lines[-mark]
	
	mark <- grep("[\\?|\\+|\\-]", lines, extended = TRUE)
	if(length(mark) > 0) lines <- lines[-mark]

	lines <- lines[grep("[^0123456789]",  lines, extended = TRUE)]

	if(datatype == 1) { # covariance matrix
		lines <- strsplit(lines, split = " ", extended = FALSE)
		lines <- lapply(lines, FUN = function(x) as.numeric(x[x!=""]))
		out <- matrix(0, nrow = length(lines), 
				 ncol = length(lines[[length(lines)]]))
		for(i in 1:nrow(out)) out[i,1:length(lines[[i]])] <- lines[[i]]
		out <- out + t(out)
		diag(out) <- diag(out) / 2
		rownames(out) <- colnames(out) <- varnames
		covlist <- list(cov = out, n.obs = n.obs)
		return(covlist)
	}
	else if(datatype == 2) { # numeric data
		out <- matrix(unlist(sapply(strsplit(lines, split = " "), 
					FUN = function(x) {
						as.numeric(x[x!=""])
					})), nrow = length(lines), byrow = TRUE)
		out <- as.data.frame(out)
		colnames(out) <- varnames
		return(out)
	}
	else if(datatype == 4) { # ordinal data
		out <- matrix(unlist(sapply(strsplit(lines, split = " "), 
					FUN = function(x) {
						as.numeric(x[x!=""])
					})), nrow = length(lines), byrow = TRUE)	
		out <- as.data.frame(out)
		for(i in 1:ncol(out)) out[,i] <- factor(out[,i], ordered = TRUE)
		colnames(out) <- varnames
		return(out)
	}
	else {
		stop("datatype not recognized, should be 1,2, or 4")
	}
}
