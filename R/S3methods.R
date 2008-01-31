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


## This file defines S3 methods

## NOTE: This file is meant to be read with 90 columns with 8 space tabs

pairs.FA <- pairs.FA.general <- pairs.FA.2ndorder <-
function(x, ...) {
	if(!is(x, "FA")) {
		stop("x must be an object of class FA or inherit from class FA")
	}
	par(pty = "s")
	Upsilon <- coef(x, "RS")
	if(ncol(Upsilon) < 2) {
		stop("pairs is sensible only for models with multiple factors")
	}
	rows <- nrow(Upsilon)
	D <- x@correlations[,,"PR"]
	Upsilon <- rbind(Upsilon, D)
	D <- diag(D)
	Phi <- x@correlations[,,"PF"]
	the.range <- c(min(c(Phi, Upsilon)), max(Upsilon))
	FC <- coef(x, "FC")
	sizes <- format(colSums(FC) / sum(FC), digits = 3)
	check.fun <- function(z, current) {
			isTRUE(all.equal(z,current, check.attributes = FALSE))
		}

	UP <- function(x, y, ...) {
		abline(v = 0, col = "gray", lty = "dashed")
		abline(h = 0, col = "gray", lty = "dashed")
		xcol <- which(apply(Upsilon, 2, check.fun, current = x))
		ycol <- which(apply(Upsilon, 2, check.fun, current = y))
		abline(v = D[xcol], col = "red", lty = "dotted")
		abline(h = D[ycol], col = "red", lty = "dotted")
		if(ncol(Upsilon) > 2) {
			other <- apply(!Upsilon[1:rows, -c(xcol, ycol), 
					drop = FALSE], 1, any)
		}
		else other <- rep(FALSE, rows)
		points(x[1:rows], y[1:rows], col = 1 + other, 
			pch = 21 - other, ...)
		}

	LP <- function(x, y, ...) {
		xcol <- which(apply(Upsilon, 2, check.fun, current = x))
		ycol <- which(apply(Upsilon, 2, check.fun, current = y))
		xpoint <- Phi[xcol,ycol]
		ypoint <- sin(acos(xpoint))
		abline(h = 0, col = "gray", lty = "dashed")
		segments(x0 = 0, y0 = 0, x1 = xpoint, y1 = 0)
		segments(x0 = 0, y0 = 0, x1 = xpoint, y1 = ypoint)
		}

  	pairs.default(Upsilon, xlim = the.range, ylim = the.range,
			labels = paste("Factor", 1:ncol(Upsilon), "\n", sizes),
			upper.panel = UP, lower.panel = LP)
}

fitted.FA <- fitted.FA.general <- fitted.FA.2ndorder <- 
function(object, ...) {
	if(!is(object, "FA")) {
		stop("object must be an object of class FA or inherit from class FA")
	}
	S <- object@manifest$cor
	if(object@restrictions@factors[1] == 1) C <- tcrossprod(coef(object))
	else C <- coef(object) %*% object@correlations[,,"PF"] %*% t(coef(object))
	return(C)
}

residuals.FA <- residuals.FA.general <- residuals.FA.2ndorder <- 
function(object, ...) {
	if(!is(object, "FA")) {
		stop("object must be an object of class FA or inherit from class FA")
	}
	S <- object@manifest$cor
	C <- fitted(object)
	R <- S - C
	return(R)
}

rstandard.FA <- rstandard.FA.general <- rstandard.FA.2ndorder <-
function(model, ...) {
	if(!is(model, "FA")) {
		stop("model must be an object of class FA or inherit from class FA")
	}
	R <- cov2cor(residuals(model))
	return(R)
}

weights.FA <- weights.FA.general <- weights.FA.2ndorder <- 
function(object, ...) {
	if(!is(object, "FA")) {
		stop("object must be an object of class FA or inherit from class FA")
	}
	if(object@method == "MLE") {
		w <- 1/tcrossprod(object@uniquenesses)
		w <- w / mean(w)
	}
	else {
		C <- fitted(object)
		communalities <- diag(C)
		w <- 1 - sweep(sweep(C^2, 1, communalities, FUN = "/"),
                                  	  2, communalities, FUN = "/")
	}
	return(w)
}

influence.FA <- influence.FA.general <- influence.FA.2ndorder <- 
function(model, ...) {
	if(!is(model, "FA")) {
		stop("model must be an object of class FA or inherit from class FA")
	}
	return(residuals(model) * weights(model))
}

predict.FA <- predict.FA.general <- predict.FA.2ndorder <- 
function(object, ...) {
	if(!is(object, "FA")) {
		stop("object must be an object of class FA or inherit from class FA")
	}
	if(is.na(object@scores[1,1])) {
		stop("cannot predict outcomes because factor scores were not generated")
	}
	return(object@scores %*% coef(object))
}

df.residual.FA <- df.residual.FA.general <- df.residual.FA.2ndorder <-
function(object, ...) {
	if(!is(object, "FA")) {
		stop("object must be an object of class FA or inherit from class FA")
	}
	return(object@restrictions@dof)
}

deviance.FA <- deviance.FA.general <- deviance.FA.2ndorder <-
function(object, ...) {
	if(!is(object, "FA")) {
		stop("object must be an object of class FA or inherit from class FA")
	}
	if(object@method == "MLE") {
		S <- model.matrix(object)
		n.obs <- object@manifest$n.obs
		llik_saturated <- -0.5 * (n.obs - 1) * (log(det(S)) + nrow(S))
		out <- llik_saturated - logLik(object)[[1]]
	}
	else    out <- -object@rgenoud$value[length(object@rgenoud$value)]
	return(out)
}

model.matrix.FA <- model.matrix.FA.general <- model.matrix.FA.2ndorder <-
function(object, ...) {
	if(!is(object, "FA")) {
		stop("object must be an object of class FA or inherit from class FA")
	}
	return(object@manifest$cor)
}
