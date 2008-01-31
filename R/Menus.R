## This file contains GUI-related functions

## NOTE: This file is meant to be read with 90 columns and 8 space tabs

FAiR_get_answer <- # everything below goes through this function
function(text, yesno = FALSE, radio_items = NULL, select = 1,
	seq_args = NULL, check_args = NULL) {
	the.env <- new.env()
	group <- ggroup(horizontal = FALSE)
	glabel(text, container = group)
	if(!is.null(radio_items) | yesno) {
		if(yesno) {
			gseparator(horizontal = FALSE, container = group)
			gobject <- gradio(c("Yes", "No"), selected = select,
					horizontal = TRUE, container = group)
		}
		else gobject <- gradio(radio_items, selected = select, container = group)
	}
	else if(!is.null(seq_args)) {
		gobject <- gspinbutton(from = seq_args$from, to = seq_args$to, 
			by = seq_args$by, value = seq_args$value, container = group)
	}
	else if(!is.null(check_args)) {
		gseparator(horizontal = FALSE, container = group)
		gobject <- gcheckboxgroup(items = check_args$items, 
					checked = check_args$checked,
					container = group)
	}
	gbasicdialog(title = if(!is.null(check_args)) 
			"Please select (multiple selections possible)" else
			"Please select one",
			widget = group,
			handler = function(h,...) {
					assign("out", value = svalue(gobject),
						envir = h$action)
			}, action = the.env,
			container = group, toolkit = guiToolkit("RGtk2"))
	if(exists("out", envir = the.env)) {
		out <- get("out", envir = the.env)
		return(out)
	}
	else stop("Cancelled", call. = FALSE)
}

FAiR_get_model <-
function() {
	items  <- c("Semi-Exploratory Factor Analysis", "Exporatory Factor Analysis",
			"Confirmatory Factor Analysis")
	answer <- FAiR_get_answer(text = "Which model would you like to estimate?",
				radio_items = items)
	if(answer == items[1]) answer <- "SEFA"
	else if(answer == items[2]) answer <- "EFA"
	else answer <- "CFA"
	return(answer)
}

FAiR_get_method <-
function() {
	items  <- c("Maximum likelihood", "Yates (1987) Weighted Least Squares")
	answer <- FAiR_get_answer(text = "How would you like to estimate the model?",
				radio_items = items)
	if(answer == items[1]) answer <- "MLE"
	else answer <- "YWLS"
	return(answer)
}

FAiR_get_algorithm <-
function() {
	text  <- "What algorithm would you like to use?"
	items <- c("call factanal() [fast]", paste("use the factanal() objective function",
		   "with RGENOUD [slower but more reliable]"), paste("CFA with zeros in",
		   "upper triangle of coefficient matrix [slower and less reliable]"))
	answer <- FAiR_get_answer(text, radio_items = items)
	return(which(items == answer))
}

FAiR_yesno <-
function(question, selected = 1) {
	answer <- FAiR_get_answer(text = question, yesno = TRUE, select = selected)
	return(answer == "Yes")
}

FAiR_get_number <-
function(text, from, to, by, value) {
	return(FAiR_get_answer(text, seq_args = list(from = from, to = to,
							by = by, value = value)))
}

FAiR_get_fix_coefficients_args <- 
function(factors, level, fixed){
	text <- paste("Would you like to put extra stipulations on the zeros at level",
			level, "?")
	items <- c("No, I just want plain zeros",
		"Encourage cohyperplanarity (revisionist Yates 1987)",
		paste(factors, "outcomes of complexity", factors - 1, 
			"for each factor (weaker Thurstone)"),
		"Unit complexity basis (revisionist Butler 1969)", 
		"Set the maximum complexity for each outcome")
	names(items) <- c("none", "quasi_Yates", "weak_Thurstone", 
			"Butler", "row_complexity")

	answer <- FAiR_get_answer(text = text, radio_items = items)
	arg_list <- formals(FAiR_fix_coefficients)
	arg_list[[names(arg_list) %in% names(which(items == answer))]] <- TRUE
	if(which(items == answer) == 5) { # row_complexity
		text <- "Would you like to set the same complexity for *all* outcomes?"
		allsame <- FAiR_yesno(question = text)
		if(allsame) {
			text <- paste("Choose the maximum complexity for all outcomes\n",
				      factors - 1, "= Thurstonian Simple Structure ...", 
				      "1 = Perfect Clustering")
			arg_list$row_complexity <- FAiR_get_number(text, from = 1, 
								by = 1, to = factors - 1, 
								     value = factors - 1)
		}
		else {
			text <- paste("What complexity would you like to use for *most*",
				      "outomes?\nYou will have the opportunity to change",
				      "this number for specific outcomes momentarily.")
			row_complexity <- FAiR_get_number(text, from = 1, by = 1,
							  to = factors, value = factors - 1)
			row_complexity <- matrix(row_complexity, nrow = nrow(fixed))
			rownames(row_complexity) <- rownames(fixed)
			text <- paste("After pressing OK, please edit the complexity for",
				      "each outome.\nWhen finished click the x at the",
				      "upper right of the editor")
			gmessage(text)
			title <- "Complexity for each test"
			row_complexity <- edit(row_complexity, title = title, 
					      edit.row.names = TRUE)
			invalid <- any(row_complexity > factors) | any(row_complexity < 1)
			while(invalid) {
				gmessage(paste("Each complexity must be between 1 and",
						factors, "inclusive. Try again."))
				row_complexity <- edit(row_complexity, title = title, 
						      edit.row.names = TRUE)
				invalid <-  any(row_complexity > factors) | 
					    any(row_complexity < 1)
			}
			arg_list$row_complexity <- c(row_complexity)
		}
	}

	text <- paste("Would you like to (default) require", factors,
			"exact zero coefficients for\neach factor at level",
			level, ", inclusive of any zeros at fixed positions?")
	zeros <- FAiR_yesno(text)
	if(!zeros) {
		zeros <- matrix(factors, ncol = factors)
		colnames(zeros) <- paste(level, "_Factor_",
					1:factors, sep = "")
		text  <- paste("After pressing OK, please edit",
			"the number of zeros required for",
			"each factor at level", level, ". Note that",
			"the number of zeros for each factor must be at least", 
			factors - 1, ", inclusive of any zeros at fixed positions.",
			"When finished, press the x in the top right of the editor.")
		invalid <- TRUE
		while(invalid) {
			gmessage(text)
			zeros <- edit(zeros, edit.row.names = FALSE)
			if(all(zeros[1,] >= (factors - 1))) invalid <- FALSE
			else {
				gmessage(paste("The number of zeros for each factor",
					"must be at least", factors - 1, ". Try again."))
				zeros[1,] <- factors
			}
		}
	}
	else zeros <- rep(factors, factors)
	arg_list$zeros <- c(zeros)
	return(as.list(arg_list))
}

FAiR_bounds_cormat <-
function(factors, level) {
	if(factors == 1) return(NULL)
	row_names <- NULL
	for(i in 1:(factors - 1)) for(j in (i+1):factors) { 
		row_names <- c( row_names, if(level == 1) 
				paste("Phi_", i, j, sep = "") else
				paste( "Xi_", i, j, sep = "") )
	}
	text <- paste("Would you like to place any non-trivial",
			"bounds on\nthe correlations among the",
			factors, "factors", "at level", level, "?\n",
			"If not, the valid interval will be (-1,1).")
	bounds <- FAiR_yesno(text, select = 2)
	if(bounds) {
		lowers <- matrix(NA_real_, nrow = factors, ncol = factors)
		lowers[upper.tri(lowers)] <- -1.0 + .Machine$double.eps
		diag(lowers) <- 1.0
		colnames(lowers) <- rownames(lowers) <- paste(
			"Factor_", 1:factors, sep = "")
		invalid <- TRUE
		while(invalid) {
			text <- paste("After pressing OK, please edit",
				"the lower bounds on the correlations",
				"among the factors at level", level, "above the",
				"diagonal and then press the x at the top right of the editor")
				
			gmessage(text)
			text <- paste("Lower bounds among factors at level", level)
			lowers <- edit(lowers, title = text, edit.row.names = TRUE)
			if(all(lowers[upper.tri(lowers)] > -1)) break
			text <- paste("All bounds on correlations must be greater",
					"than -1.0. Try again.")
			gmessage(text)
		}

		uppers <- lowers
		uppers[upper.tri(uppers)] <- 1.0 - .Machine$double.eps
		while(invalid) {
			text <- paste("After pressing OK, please edit",
				"the upper bounds on the correlations",
				"among the factors at level", level, "above the",
				"diagonal and then press the x at the top right of the editor")
			gmessage(text)
			text <- paste("Upper bounds among factors at level", level)
			uppers <- edit(uppers, title = text, edit.row.names = TRUE)
			if(all(uppers[upper.tri(uppers)] < 1)) break
			text <- paste("All bounds on correlations must be less",
					"than 1.0. Try again.")
			gmessage(text)
		}
		Domains <- cbind(lowers[upper.tri(lowers)], uppers[upper.tri(uppers)])
		rownames(Domains) <- row_names
	}
	else if(factors == 2) { 
		Domains <- cbind(-1 + .Machine$double.eps, 1 - .Machine$double.eps)
		rownames(Domains) <- if(level == 1) "Phi_12" else "Xi_12"
	}
	else {
		Domains <- cbind(-1 + .Machine$double.eps, rep(1 - .Machine$double.eps,
				0.5 * factors * (factors - 1)))
		rownames(Domains) <- row_names
	}
	return(Domains)
}

FAiR_bounds_coefficients <-
function(fixed, level) {
	thing <- if(level == 1) 1 else if(ncol(fixed) == 1) 2 else 3
	text <- paste("Would you like to place any non-trivial bounds\non",
			"any of the (free) coefficients at level", level, "?\n",
			"If not, the valid inteval will be", 
			switch(thing, "[-1.5, 1.5]", "[-1.0, 1.0]", "[-1.5, 1.5]"), ".")
	bounds <- FAiR_yesno(text, select = 2)
	if(bounds) {
		text <- paste("What number would you like to use as a lower bound",
				"for *most* of the coefficients?\nYou will have an",
				"opportunity to alter this number for specific",
				"coefficients momentarily.")
		num <- FAiR_get_number(text, from = switch(thing, -10, -1.0, -10),
					to = switch(thing, 1.5, 1.0, 1.5), by = .1,
					value = switch(thing, -1.5, -1.0, -1.5))
		lowers <- matrix(num, nrow = nrow(fixed), ncol = ncol(fixed))
		rownames(lowers) <- rownames(fixed)
		colnames(lowers) <- colnames(fixed)
		lowers[!is.na(fixed)] <- fixed[!is.na(fixed)]
		if(level == 1) {
			text <- paste("After pressing OK, please edit the lower bounds",
			  		"on the coefficients at level 1.", 
					"Note that it is theoretically",
					if(ncol(fixed) > 1) "possible" else "impossible",
					"for a coefficient to be less than -1.0.\n",
					"When finished press the x at the top right of the editor.")
		}
		else {
			text <- paste("After pressing OK, please edit the lower bounds",
			  		"on the coefficients at level 2.", 
					"Note that it is theoretically",
					if(ncol(fixed) > 1) "possible" else "impossible",
					"for a coefficient to be less than -1.0.\n",
					"When finished press the x at the top right of the editor.")
		}
		invalid <- TRUE
		while(invalid) {
			gmessage(text)
			lowers <- edit(lowers, edit.row.names = TRUE)
			if( (ncol(fixed) == 1) && any(lowers < -1) ) {
				gmessage("All bounds must be greater than -1.0")
			}
			else break
		}

		text <- paste("What number would you like to use as an upper bound",
				"for *most* of the coefficients?\nYou will have an",
				"opportunity to alter this number for specific",
				"coefficients momentarily.")
		num <- FAiR_get_number(text, from = switch(thing, -1.5, -1.0, -1.5), 
					to = switch(thing, 10, 1.0, 10), by = .1,
					value = switch(thing, 1.5, 1.0, 1.5))
		uppers <- matrix(num, nrow = nrow(fixed), ncol = ncol(fixed))
		uppers[!is.na(fixed)] <- fixed[!is.na(fixed)]
		rownames(uppers) <- rownames(fixed)
		colnames(uppers) <- colnames(fixed)
		if(level == 1) {
			text <- paste("After pressing OK, please edit the upper bounds",
			  		"on the coefficients at level 1.", 
					"Note that it is theoretically",
					if(ncol(fixed) > 1) "possible" else "impossible",
					"for a coefficient to be greater than 1.0.\n",
					"When finished press the x at the top right of the editor.")
		}
		else {
			text <- paste("After pressing OK, please edit the upper bounds",
			  		"on the coefficients at level 2.", 
					"Note that it is theoretically",
					if(ncol(fixed) > 1) "possible" else "impossible",
					"for a coefficient to be greater than 1.0.\n",
					"When finished press the x at the top right of the editor.")
		}
		while(invalid) {
			gmessage(text)
			uppers <- edit(uppers, edit.row.names = TRUE)
			if( (ncol(fixed) == 1) && any(lowers > 1) ) {
				gmessage("All bounds must be less than 1.0")
			}
			else break
		}
		Domains <- cbind(c(lowers), c(uppers))
	}
	else {
		bound <- switch(thing, 1.5, 1.0, 1.5)
		Domains <- cbind(-bound, rep(bound, prod(dim(fixed))))
	}
	row_names <- NULL
	for(i in 1:ncol(fixed)) row_names <- c(row_names, 
						paste(rownames(fixed), "_", i, sep = ""))
	rownames(Domains) <- row_names
	return(Domains)
}

FAiR_peg_coefficients <-
function(coefs, level) {
	text <- paste("After pressing OK, please edit the following matrix",
			"to restrict coefficients to specified values at level",
			level, ".\nLeave any *unrestricted* coefficients as NA.\n",
			"When finished press the x at the top right of the editor.")
	gmessage(text)
	text <- paste("Matrix of restricted coefficients")
	rows <- nrow(coefs)
	cols <- ncol(coefs)
	coefs <- edit(coefs, title = text, edit.row.names = TRUE)
	invalid <- cols == 1
	while(invalid) {
		if(all(coefs >= -1.0, na.rm = TRUE) & 
		   all(coefs <=  1.0, na.rm = TRUE)) break
		gmessage("All coefficients must be on the [-1.0,1.0] interval, try again")
		coefs <- edit(coefs, title = text, edit.row.names = TRUE)
	}
	return(coefs[1:rows,1:cols,drop = FALSE])
}

FAiR_criterionator_extraction <-
function(levels, factors, method, criteria) {
	if(!is.null(criteria)) {
		if(!is.list(criteria)) {
			stop("'criteria' must be a list or left unspecified")
		}
		else if(length(criteria) > 0) {
			for(i in 1:length(criteria)) {
				if(is.character(criteria[[i]])) { 
					criteria[[i]] <- get(paste("FAiR_criterion_",
								criteria[[i]], sep = ""))
				}
				else if(!is.function(criteria[[i]])) {
					stop("criteria must be a list of character",
					     " strings or functions, it is usually best",
					     " to leave criteria unspecified")
				}
			}
		}
		else if(method == "MLE")  {
			criteria[[1]] <- FAiR_criterion_llik
			names(criteria) <- "llik"
		}
		else if(method == "YWLS") {
			criteria[[1]] <- FAiR_criterion_YWLS
			names(criteria) <- "YWLS"
		}
		return(criteria)
	}
	items <-  c(
		paste("Reference factors at level 2 must have more effective",
		      "variance than do the primary factors at level 1"),
		paste("Primary factors at level 2 must have more effective",
		      "variance than do the primary factors at level 1"),
		paste("Reference factors at level 1 must have more effective",
		      "variane than does the battery as a whole"),
		paste("Primary factors at level 1 must have more effective",
		      "variance than does the battery as a whole"),
		"No suppressor variables at level 2",
		"No suppressor variables at level 1",
		paste("Reference factors at level 2 must have more generalized",
		      "variance than do the primary factors at level 2"),
		paste("Reference factors at level 1 must have more generalized",
		      "variance than do the primary factors at level 1"),
		paste("Tests in hyperplanes have more effective variance than",
		      "does the battery as a whole"))
	names(items) <- c("evRF_2nd", "evPF_2nd", "evRF_1st", "evPF_1st", 
			  "no_suppressors_2nd", "no_suppressors_1st", 
			  "dets_2nd", "dets_1st", "cohyperplanarity")
	if( (levels == 1) | (factors[2] == 1) ) {
		mark <- c(3, 4, 6, 8, 9)
		items <- items[mark]
		checked <- c(1, 0, 1, 0, 0)
	}
	else {
		checked <- c(0, 0, 1, 0, 1, 1, 0, 0, 0)
	}
	checked <- as.logical(checked)

	if(method == "MLE") {
		text <-  paste( "Please indicate which, if any, lexical criteria you",
				"would like to use.\nThe log-likelihood will",
				"automatically be added as the last, or only, criterion.")
	}
	else if(method == "YWLS") {
		text <-  paste( "Please indicate which, if any, lexical criteria you",
				"would like to use.\nThe WLS criterion in Yates (1987)",
				"will be added as the last, or only, criterion")
	}
	else stop("This should not have happened. Please notify maintainer.")

	answer <- FAiR_get_answer(text, check_args = list(items = items, 
				  checked = checked))
	if(length(answer) > 0) {
		answer_names <- names(answer)
		criteria <- list()
		threshold <- NA_real_
		for(i in 1:length(answer_names)) {
			mesh <- paste("FAiR_criterion_", answer_names[i], sep = "")
			criteria[[i]] <- get(mesh)
			names(criteria)[[i]] <- answer_names[i]
			if(answer_names[i] %in% c("no_suppressors_2nd", 
						  "no_suppressors_1st")) {
				if(is.na(threshold)) {
					text <- paste("Select the minimum acceptable",
						"factor contribution coefficient\nsuch",
						"that more negative values imply a",
						"suppressor variable.")
					threshold <- FAiR_get_number(text, from = -1, 
							to = 0, by = .001, value = -.01)
				}
				formals(criteria[[i]])$threshold <- threshold
			}
		}
	}
	else criteria <- list()
	criteria[[length(criteria) + 1]] <- if(method == "MLE") FAiR_criterion_llik else
                                                                FAiR_criterion_YWLS
	names(criteria)[length(criteria)] <- if(method == "MLE") "llik" else "YWLS"
	return(criteria)
}

make_restrictions <- 
function(factors, model, method, fixed, covmat, criteria = NULL) {
	require(gWidgetsRGtk2)
	options(guiToolkit = "RGtk2")
	flush.console()
	cat("Note: The GUI window may have popped up behind your other windows.\n")
	flush.console()
	if(missing(model))  model  <- FAiR_get_model()
	else                model  <- match.arg(model,  eval(formals(Factanal)$model))

	if(missing(method)) method <- FAiR_get_method()
	else                method <- match.arg(method, eval(formals(Factanal)$method))

	if(missing(covmat)) stop("covmat must be specified")
	S <- FAiR_parse(covmat = covmat, robust = FALSE, seeds = 12345)$cor

	if(model == "EFA") {
		levels <- 1
		if(missing(factors)) {
			text <- "How many factors should be extracted?"
			factors <- FAiR_get_number(text, from = 1, to = floor(nrow(S)/2),
							value = 3, by = 1)
		}
		else stopifnot(is.numeric(factors), length(factors) <= 2,
				factors == as.integer(factors), factors[1] > 0)
		if(length(factors) == 1) factors <- c(factors, 0)
		else factors[2] <- 0

		p <- nrow(S)
		dof <- as.integer(0.5 * ((p - factors[1])^2 - p - factors[1]))
		Phi <- diag(factors[1])
		attributes(Phi)$ev <- 1
		Theta2 <- diag(p)
		if(method == "MLE") {
			algorithm <- FAiR_get_algorithm()
			if(algorithm < 3) {
				Domains <- cbind(0, rep(1, p))
				rownames(Domains) <- rownames(S)
				restrictions <- new("restrictions.factanal", 
						factors = factors, nvars = p,
						dof = dof, Domains = Domains,
						model = model, method = method,
						fast = algorithm == 1)
			}
			else {
				top <- diag(factors)
				top[lower.tri(top, diag = TRUE)] <- NA_real_
				fixed <- rbind(top, matrix(NA_real_, nrow = p - factors, 
								     ncol = factors) )
				rownames(fixed) <- rownames(S)
				Domains <- cbind(-1, rep(1, sum(is.na(fixed))))
				Domains[1,1] <- 0
				Domains <- rbind(Domains, cbind(0, rep(1, p)))
				beta_list <- list(beta = fixed, free = c(is.na(fixed)),
						num_free = sum(is.na(fixed)))
				restrictions <- new("restrictions.orthonormal",
						factors = factors, nvars = nrow(Domains),
						dof = dof, Domains = Domains,
						model = model, method = method,
						Phi = Phi, beta = beta_list,
						Theta2 = list(Theta2 = Theta2),
						criteria = list(FAiR_criterion_llik))
			}
		}
		else { # YWLS
			top <- diag(factors)
			top[lower.tri(top, diag = TRUE)] <- NA_real_
			fixed <- rbind(top, matrix(NA_real_, nrow = p - factors, 
							     ncol = factors) )
			rownames(fixed) <- rownames(S)
			Domains <- cbind(-1, rep(1, sum(is.na(fixed))))
			Domains[1,1] <- sqrt(.Machine$double.eps)
			Domains <- rbind(Domains, cbind(0, rep(1, nrow(S))))
			beta_list <- list(beta = fixed, free = c(!is.na(fixed)))
			restrictions <- new("restrictions.orthonormal",
					factors = factors, nvars = nrow(Domains),
					dof = dof, Domains = Domains,
					model = model, method = method,
					Phi = Phi, beta = beta_list,
					Theta2 = Theta2,
					criteria = list(FAiR_criterion_llik))
		}
	}
	else { # SEFA / CFA
		if(missing(factors)) {
			text <- "How many factors should be extracted at level 1?"
			factors <- FAiR_get_number(text, from = 1, to = floor(nrow(S)/2),
							value = 2, by = 1)
		}
		else stopifnot(is.numeric(factors), all(factors == as.integer(factors)), 
				factors[1] > 0)

		if(factors[1] >= 3) {
			text <- paste("Would you like to estimate a simultaneous",
					"second-order model?")
			second <- FAiR_yesno(text)
			levels <- 1 + second

			if(levels > 1 & (length(factors) == 1)) {
				text <- "How many factors should be extracted at level 2?"
				factors <- if(factors[1] < 5)  c(factors, 1) else
						c(factors, FAiR_get_number(text, from = 1,
						to = floor(factors/2), by = 1, value = 1))
			}
			else if(second && (length(factors) == 1)) factors <- c(factors, 1)
			else if(length(factors) == 1) factors <- c(factors, 0)
		}
		else { 
			levels <- 1
			if(length(factors) == 1) factors <- c(factors, 0)
		}

		if(missing(fixed)) {
			fixed <- matrix(NA_real_, nrow = nrow(S), ncol = factors[1])
			rownames(fixed) <- rownames(S)
			colnames(fixed) <- paste("Factor_", 1:factors[1], sep = "")

			if(factors[2] > 0) {
				fixed2 <- matrix(NA_real_, factors[1], factors[2])
				rownames(fixed2) <- paste("1_Factor_", 1:factors[1], 
								sep = "")
				colnames(fixed2) <- paste("2_Factor_", 1:factors[2], 
								sep = "")
			}
		}
		else if(is.list(fixed)) {
			if(levels == 1) {
				stop("you passed 'fixed' as a list but specified a",
					"single-equation model, please respecify")
			}
			fixed2 <- fixed[[2]]
			fixed  <- fixed[[1]]

			if(nrow(fixed) != nrow(S)) {
				stop("'fixed[[1]]' must have as many rows as 'S'")
			}
			if(ncol(fixed) != factors[1]) {
				stop("'fixed[[1]]' must have as many columns as there",
					" are factors at level 1")
			}
			if(is.null(rownames(fixed))) rownames(fixed) <- rownames(S)

			if(nrow(fixed2) != factors[1]) {
				stop("'fixed[[2]]' must have as many rows as there are",
					"factors at level 1")
			}
			else if(ncol(fixed2) != factors[2]) {
				stop("'fixed[[2]]' must have as many columns as there",
					"are factors are level 2")
			}
			if(is.null(rownames(fixed2))) {
				rownames(fixed2) <- paste("1_Factor_", 1:factors[1], 
							sep = "")
			}
			if(is.null(colnames(fixed2))) {
				colnames(fixed2) <- paste("2_Factor_", 1:factors[2], 
							sep = "")
			}
		}
		else if(is.matrix(fixed)) {
			if(nrow(fixed) != nrow(S)) {
				stop("'fixed' must have as many rows as 'S'")
			}
			if(ncol(fixed) != factors[1]) {
				stop("'fixed' must have as many columns as there are",
					" factors at level 1")
			}
			if(is.null(rownames(fixed))) rownames(fixed) <- rownames(S)

			if(factors[2] > 0) {
				fixed2 <- matrix(NA_real_, factors[1], factors[2])
				rownames(fixed2) <- paste("1_Factor_", 1:factors[1], 
								sep = "")
				colnames(fixed2) <- paste("2_Factor_", 1:factors[2], 
								sep = "")
			}
		}
		else {
			stop("'fixed' must be a matrix, a list of two matrices, or ",
				"left unspecified (recommended)")
		}

		if(levels > 1) {
			if(factors[2] > 1) Domains <- FAiR_bounds_cormat(factors[2], 2)
			else Domains <- matrix(NA_real_, nrow = 0, ncol = 2)
			if(model == "SEFA") {
				text <- paste("Would you like to peg any coefficients",
						"to exact values at level 2?")
				pegs <- FAiR_yesno(text, selected = 2)
				if(pegs) fixed2 <- FAiR_peg_coefficients(fixed2, 2)
			}
			else if(all(is.na(fixed2))) {
				fixed2 <- FAiR_peg_coefficients(fixed2, 2) # CFA
			}

			temp_Domains <- FAiR_bounds_coefficients(fixed2, level = 2)
			pegged2 <- !is.na(fixed2)
			if(any(pegged2)) temp_Domains <- temp_Domains[!c(pegged2),]
			Delta_select <- rep(FALSE, nrow(Domains))
			Domains <- rbind(Domains, temp_Domains)
			Delta_select <- c(Delta_select, rep(TRUE, sum(!pegged2)))
			Delta_list <- list(Delta = fixed2, free = c(is.na(fixed2)),
						select = Delta_select) # added to later

			if( (model == "SEFA") && (factors[2] > 1) ) {
				fix_Delta_args <- FAiR_get_fix_coefficients_args(
							factors = factors[2], level = 2,
							fixed = fixed2)
				Delta_list$fix_Delta_args <- fix_Delta_args
			}
		}
		else Domains <- FAiR_bounds_cormat(factors[1], 1)

		if(model == "SEFA") {
			text <- if(levels == 1) paste("Would you like to peg any",
						"coefficients to exact values?") else
				paste("Would you like to peg any coefficients at level",
					"1 to exact values?")
			pegs <- FAiR_yesno(text, selected = 2)
			if(pegs) fixed <- FAiR_peg_coefficients(fixed, 1)
		}
		else if(all(is.na(fixed))) {
			fixed <- FAiR_peg_coefficients(fixed, 1)
		}

		temp_Domains <- FAiR_bounds_coefficients(fixed, level = 1)
		pegged <- !is.na(fixed)
		if(any(pegged)) temp_Domains <- temp_Domains[!c(pegged),]
		beta_select <- rep(FALSE, nrow(Domains))
		Domains <- rbind(Domains, temp_Domains)
		beta_select <- c(beta_select,   rep(TRUE, sum(!pegged)), 
						rep(FALSE, nrow(S)))
		if(model == "SEFA") {
			fix_beta_args <- FAiR_get_fix_coefficients_args(factors[1], 1, fixed)
			indeterminate <- FAiR_indeterminator("SEFA", factors,
							zeros1 = fix_beta_args$zeros,
							zeros2 = if(levels == 1) 0 else 
								 if(factors[2] == 1) 0 else
								fix_Delta_args$zeros,
							nonzeros1 = sum(fixed  != 0, na.rm = TRUE),
							nonzeros2 = sum(fixed2 != 0, na.rm = TRUE))
			dof <- 0.5 * nrow(S) * (nrow(S) + 1) - nrow(Domains) - nrow(S) + 
				sum(fix_beta_args$zeros) - sum(fixed == 0, na.rm = TRUE) +
				sum(Domains[,1] == Domains[,2])
			if(!is.na(fix_beta_args$row_complexity[1])) {
				dof <- dof - sum(fix_beta_args$zeros)
				if(length(fix_beta_args$row_complexity) == 1) {
					dof <- dof + nrow(S) * (factors[1] - 
								fix_beta_args$row_complexity)
				}
				else dof <- dof + sum(factors[1] - fix_beta_args$row_complexity)
					
			}
			if(levels > 0) {
				if(factors[2] > 1) {
					dof <- dof + sum(fix_Delta_args$zeros) - 
						     sum(fixed2 == 0, na.rm = TRUE)
					if(!is.na(fix_Delta_args$row_complexity[1])) {
						dof <- dof - sum(fix_Delta_args$zeros)
						if(length(fix_Delta_args$row_complexity) == 1) {
							dof <- dof + factors[1] * (factors[2] - 
								   fix_Delta_args$row_complexity)
						}
						else dof <- dof + sum(factors[2] - 
								   fix_Delta_args$row_complexity)
					}
				}
			}
			beta_list <- list(beta = fixed, free = c(!pegged), 
					select = beta_select, 
					fix_beta_args = fix_beta_args)
		}
		else { # CFA
			indeterminate <- FAiR_indeterminator("CFA", factors,
							zeros1 = sum(fixed == 0, na.rm = TRUE),
							zeros2 = sum(fixed2 == 0, na.rm = TRUE),
							nonzeros1 = sum(fixed != 0, na.rm = TRUE),
							nonzeros2 = sum(fixed2 != 0, na.rm = TRUE))
			if(indeterminate == "determined") {
				temp_fixed <- fixed
				temp_fixed[is.na(temp_fixed)] <- runif(sum(is.na(temp_fixed)), 
									max = 0.5)
				check1 <- FAiR_check_coefficients(temp_fixed)
				if(check1 != 1) indeterminate <- "Howe1"
				else if(levels > 1) {
					temp_fixed <- fixed2
					temp_fixed[is.na(temp_fixed)] <- runif(sum(is.na(temp_fixed)),
										max = 0.5)
					check2 <- FAiR_check_coefficients(temp_fixed)
					if(check2 != 1) indeterminate <- "Howe2"
				}
			}
			dof <- 0.5 * nrow(S) * (nrow(S) + 1) - nrow(Domains) - nrow(S)
			beta_list <- list(beta = fixed, free = c(!pegged), 
					select = beta_select)
		}

		dof <- as.integer(dof)
		if(dof < 0) {
			stop("negative degrees of freedom", 
			     "estimate fewer factors or impose more restrictions")
		}

		if(indeterminate != "determined") {
			if(indeterminate == "minimal") {
				text <- paste("It appears that the restrictions only minimally",
						"satisfy the theorem in Howe (1955) on rotational",
						"indeterminacy.\n A SEFA model can be estimated in",
						"this case, but it is not clear that the mode of the",
						"objective function is unique.\n It is recommended",
						"to press 'Cancel' and impose at least one more",
						"restriction on the model.\n Or you can press 'OK'",
						"to continue if you are sure the mapping rule",
						"or bounds bind.")
			}
			else if(indeterminate == "Howe1") {
				text <- paste("It appears that the restrictions do not satisfy the",
						"rank condition in Howe (1955) on rotational",
						"indeterminacy at level 1.\n It is recommended to",
						"'Cancel' and rectify this situation but you can",
						"press 'OK' to continue if you are sure that",
						"rotational indeterminacy is eliminated.")
				warning("It appears the rank condition is not satisfied at level 1.")
			}
			else if(indeterminate == "Howe2") {
				text <- paste("It appears that the restrictions do not satisfy the",
						"rank condition in Howe (1955) on rotational",
						"indeterminacy at level 2.\n It is recommended to",
						"'Cancel' and rectify this situation but you can",
						"press 'OK' to continue if you are sure that",
						"rotational indeterminacy is eliminated.")
				warning("It appears the rank condition is not satisfied at level 2.")
			}
			else {
				text <- paste("There appear to be too few restrictions",
						"to eliminate rotational indeterminacy",
						"according\nto the theorem in Howe (1955).",
						"Click OK to continue anyway or press Cancel.")
				warning("There appear to be too few restrictions to",
					" eliminate rotational indeterminacy.")
			}
			stopifnot(gbasicdialog(widget = glabel(text)))
		}
		criteria  <- FAiR_criterionator_extraction(levels, factors, 
                                                           method, criteria)
		Theta2_select <- rep(FALSE, nrow(Domains))
		Domains <- rbind(Domains, cbind(0, rep(1, nrow(S))))
		Theta2_select <- c(Theta2_select, rep(TRUE, nrow(S)))
		Theta2_list <- list(Theta2 = diag(nrow(S)), diagonal = TRUE,
					select = Theta2_select)
		diag(Theta2_list$Theta2) <- NA_real_
		colnames(Domains) <- c("lower", "upper")
		if(levels == 1) {
			restrictions <- new("restrictions.1storder",
						factors = factors, nvars = nrow(Domains),
						dof = dof, Domains = Domains,
						model = model, method = method,
						Phi = diag(rep(0.5, factors[1])),
						beta = beta_list, Theta2 = Theta2_list,
						criteria = criteria)
		}
		else if(factors[2] == 1) {
			Delta_list$select <- NULL
			Delta_list$num_free <- sum(Delta_list$free)
			restrictions <- new("restrictions.general",
						factors = factors, nvars = nrow(Domains),
						dof = dof, Domains = Domains,
						model = model, method = method,
						Phi = diag(rep(0, factors[1])),
						Delta = Delta_list,
						beta = beta_list, Theta2 = Theta2_list,
						criteria = criteria)
		}
		else {
			# need to pad Delta_list$select with FALSE
			Delta_list$select <- c(Delta_list$select, rep(FALSE, 
					       nrow(Domains) - length(Delta_list$select)))
			restrictions <- new("restrictions.2ndorder",
						factors = factors, nvars = nrow(Domains),
						dof = dof, Domains = Domains,
						model = model, method = method,
						Xi = diag(rep(0.5, factors[2])),
						Delta = Delta_list,
						Phi = diag(rep(0, factors[1])),
						beta = beta_list, Theta2 = Theta2_list,
						criteria = criteria)
		}
	}
	return(restrictions)
}

make_criteria <-
function(factors, weights, FAobject) {
	require(gWidgetsRGtk2)
	options(guiToolkit = "RGtk2")
	cat("Note: The GUI window may have popped up behind your other windows.\n")
	flush.console()

	criteria <- list()

	# No factor collapse criterion
	text <- paste("Select the minimum acceptable effective variance for the",
			"correlation matrix among factors.\nSetting this threshold",
			"too close to zero risks factor collapse.")
	threshold     <- FAiR_get_number(text, from = 0, to = .9, by = .01, value = .25)
	criteria[[1]] <- FAiR_criterion_no_factor_collapse
	formals(criteria[[1]])$threshold <- threshold

	# Other restriction criteria
	items <- c("Limit primary factor correlations", 
		paste("Reference factors must have more effective",
		      "variance than does the battery as a whole"),
		paste("Primary factors must have more effective",
		      "variance than does the battery as a whole"),
		"No suppressor variables", 
		paste("Reference factors must have more generalized",
		      "variance than do the primary factors"),
		"Positive manifold")
	names(items) <- paste("FAiR_criterion_", c("limit_correlations", "ev_RF", "ev_PF",
			      "dets_1st", "no_suppressors", "positive_manifold"), sep = "")
				
	text <-  paste( "Please indicate which, if any, additional lexical criteria you",
			"would like to use.\nYou will have an opportunity to select",
			"the ultimate criterion momentarily")
	answer <- FAiR_get_answer(text, check_args = list(items = items, 
					checked = c(FALSE, FALSE, TRUE, FALSE)))
	if(length(answer) > 0) {
		answer_names <- names(answer)
		if(names(items)[1] %in% answer_names) { # limit inter-factor correlations
			criteria[[length(criteria)+1]] <- get(names(items)[1])

			# set lower bound
			text  <- paste("Select the minimum acceptable correlation",
					"among primary factors")
			lower <- FAiR_get_number(text, from = -1, to = 1, 
						 by = .01, value = -1)
			formals(criteria[[length(criteria)]])$lower <- lower

			# set upper bound
			text  <- paste("Select the maximum acceptable correlation",
					"among primary factors")
			upper <- FAiR_get_number(text, from = lower, to = 1, 
						  by = .01, value = 1)
			formals(criteria[[length(criteria)]])$upper <- upper

			if(lower == 0 && upper == 0) {
				warning("Rotate() does not do orthogonal rotation.",
					"Ignoring bounds on inter-factor correlations.")
			}
		}
		if(names(items)[2] %in% answer_names) { # ev_RF
			criteria[[length(criteria)+1]] <- FAiR_criterion_dets
			C <- fitted(FAobject) + diag(FAobject@uniquenesses)
			formals(criteria[[length(criteria)]])$threshold <- det(C)^(1/ncol(C))				
		}
		if(names(items)[3] %in% answer_names) { # ev_PF
			criteria[[length(criteria)+1]] <- FAiR_criterion_dets
			C <- fitted(FAobject) + diag(FAobject@uniquenesses)
			formals(criteria[[length(criteria)]])$threshold <- det(C)^(1/ncol(C))
		}
		if(names(items)[4] %in% answer_names) { # no suppressors
			criteria[[length(criteria)+1]] <- FAiR_criterion_no_suppressors
			text <- paste("Select the minimum acceptable factor contribution",
				"coefficient\nsuch that lower values indicate a",
				"suppressor variable.")
			threshold <- FAiR_get_number(text, from = -1, to = 0,
					by = .001, value = -.01)
			formals(criteria[[length(criteria)]])$threshold <- threshold
		}
		if(names(items)[5] %in% answer_names) { # dets
			criteria[[length(criteria)+1]] <- FAiR_criterion_dets
		}
		if(names(items)[6] %in% answer_names) { # positive manifold
			criteria[[length(criteria)+1]] <- FAiR_criterion_positive_manifold
			text <- paste("Select the minimum acceptable reference structure",
					"correlation.\nA slightly negative number is",
					"recommended rather than zero")
			threshold <- FAiR_get_number(text, from = -1, to = 0, 
							by = .01, value = -.1)
			formals(criteria[[length(criteria)]])$threshold <- threshold
		}
	}

	# Ultimate criterion
	text   <- paste("Please indicate which analytic criterion should be the ultimate",
			"criterion during the lexical minimization.")
	items  <- c(paste("minimaximin: The maximum of the minimum squared reference",
			"structure correlations for each outcome."), paste("phi:", 
			"the criterion proposed by Thurstone"), paste("varphi:",
			"a generalization of phi using weighted sums"),
			"LS: Loading Simplicity Index proposed by Lorenzo-Seva")
	names(items) <- c("minimaximin", "phi", "varphi", "LS")
	answer <- FAiR_get_answer(text, radio_items = items, select = 3)
	if(answer == items["phi"]) {
		criteria[[length(criteria)+1]] <- FAiR_criterion_phi
		text <- paste("In calculating phi, the reference structure correlations",
				"will be raised to the power 2/c.\nPlease choose c,",
				"which need not be an integer but is usually taken to be",
				"1.0 by Thurstone.")
		c <- FAiR_get_number(text, from = 1, to = 20, by = .1, value = 1)
		formals(criteria[[length(criteria)]])$c <- c
	}
	else if(answer == items["minimaximin"]) {
		criteria[[length(criteria)+1]] <- FAiR_criterion_minimaximin
	}
	else if(answer == items["LS"]) {
		criteria[[length(criteria)+1]] <- FAiR_criterion_LS
	}
	else {
		criteria[[length(criteria)+1]] <- FAiR_criterion_varphi
		if(!is.null(weights)) {
			formals(criteria[[length(criteria)]])$weights <- weights
			return(criteria)
		}
		text  <- paste("What kind of weights would you like to use?")
		items <- c("Dynamic weights", "User-specified weights")
		weights <- FAiR_get_answer(text, radio_items = items)
		if(weights == items[[2]]) {
			weight  <- 1
			weights <- rep(0, factors - 1)
			for(i in 1:length(weights)) {
				weight  <- FAiR_get_number(paste("Select weight", i), 
					from = 0, to = weight, by = .01, value = weight)
				weights[[i]] <- weight
			}
			formals(criteria[[length(criteria)]])$weights <- weights
		}
	}
	return(criteria)
}
