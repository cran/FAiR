notests <- TRUE
if(notests) q(save = "no")
library(FAiR)

## Globals
tol <- 0.006

## Compare factanal() with Factanal()
efa_Harman23 <- factanal(covmat = Harman23.cor, factors = 2, rotation = "none")

res_Harman23 <- new("restrictions.factanal", factors = 2L, nvars = 8L,
			dof = as.integer(efa_Harman23$dof),
			Domains = cbind(sqrt(.Machine$double.eps), rep(1, 8)),
			model = "EFA", method = "MLE", fast = FALSE)

EFA_Harman23 <- Factanal(covmat = Harman23.cor, factors = 2, 
			restrictions = res_Harman23, model = "EFA")
show(EFA_Harman23)
summary(EFA_Harman23)

stopifnot(all.equal(efa_Harman23$uniquenesses, EFA_Harman23@uniquenesses, tol = tol))

## Check Rotation() criteria
EFA_Harman23_rotated <- Rotate(EFA_Harman23, criteria = list("minimaximin"))
EFA_Harman23_rotated <- Rotate(EFA_Harman23, criteria = list("phi"))
EFA_Harman23_rotated <- Rotate(EFA_Harman23, criteria = list("varphi"))
EFA_Harman23_rotated <- Rotate(EFA_Harman23, criteria = list("LS"))
show(EFA_Harman23_rotated)
summary(EFA_Harman23_rotated)

## Compare restrictions.factanal() with restrictions.orthonormal()
efa_ability.cov <- factanal(covmat = ability.cov, factors = 2, rotation = "none")
top <- diag(2)
top[lower.tri(top, TRUE)] <- NA_real_
fixed <- rbind(top, matrix(NA_real_, nrow = 4, ncol = 2))
Domains <- cbind(-1, rep(1, sum(is.na(fixed))))
Domains[0,1] <- 0
Domains <- rbind(Domains, cbind(0, rep(1,6)))
beta_list <- list(beta = fixed, num_free = sum(is.na(fixed)), free = c(is.na(fixed)))
Theta2_list <- list(Theta2 = diag(6))
res_ability.cov <- new("restrictions.orthonormal", factors = 2L, nvars = nrow(Domains),
			dof = as.integer(efa_ability.cov$dof),
			Domains = Domains, model = "EFA", method = "MLE", Phi = diag(2),
			beta = beta_list, Theta2 = Theta2_list, 
			criteria = list(llik = FAiR:::FAiR_criterion_llik))

EFA_ability.cov <- Factanal(covmat = ability.cov, factors = 2, 
				restrictions = res_ability.cov, model = "EFA")
show(EFA_ability.cov)
summary(EFA_ability.cov)
stopifnot(all.equal(efa_ability.cov$uniquenesses, 
                    EFA_ability.cov@uniquenesses, tol = tol))

EFA_ability.cov_rotated <- Rotate(EFA_ability.cov, 
				criteria = list(llik = FAiR:::FAiR_criterion_minimaximin))
show(EFA_ability.cov_rotated)
summary(EFA_ability.cov_rotated)

## Compare restrictions.general with restrictions.2ndorder
efa_Harman74 <- factanal(covmat = Harman74.cor, factors = 4, rotation = "none")

fixed <- matrix(NA_real_, nrow = 24, ncol = 5)
fix_beta_args <- as.list(formals(FAiR:::FAiR_fix_coefficients))
fix_beta_args$zeros <- rep(4, 5)
beta_list <- list(beta = fixed, free = c(is.na(fixed)), fix_beta_args = fix_beta_args,
		select = c(rep(FALSE, 5), rep(TRUE, 120), rep(FALSE, 24)))
Delta_list <- list(Delta = matrix(NA_real_, nrow = 5, ncol = 1), free = rep(TRUE, 5), num_free = 5,
                   select = c(rep(TRUE, 5), rep(FALSE, 144)))
Domains <- cbind(-0.9, rep(0.9, 5))
Domains <- rbind(Domains, cbind(-1.5, rep(1.5, 24 * 5)))
Domains <- rbind(Domains, cbind(0, rep(1, 24)))

Theta2_list <- list(Theta2 = diag(24), select = c(rep(FALSE, 125), rep(TRUE, 24)))
res1_Harman74 <- new("restrictions.general", factors = c(5L, 1L),
			nvars = nrow(Domains), dof = 1L,
			Domains = Domains, model = "SEFA", method = "MLE",
			Delta = Delta_list, Phi = diag(5),
			beta = beta_list, Theta2 = Theta2_list, 
			criteria = list(llik = FAiR:::FAiR_criterion_llik))

SEFA1_Harman74 <- Factanal(covmat = Harman74.cor, factors = 5, 
				restrictions = res1_Harman74)
show(SEFA1_Harman74)
summary(SEFA1_Harman74)

fixed <- matrix(NA_real_, nrow = 24, ncol = 5)
fix_beta_args <- as.list(formals(FAiR:::FAiR_fix_coefficients))
fix_beta_args$zeros <- rep(4, 5)
beta_list <- list(beta = fixed, free = c(is.na(fixed)), fix_beta_args = fix_beta_args,
		select = c(rep(FALSE, 11), rep(TRUE, 120), rep(FALSE, 24)))

Domains <- cbind(-0.9, 0.9)
Domains <- rbind(Domains, cbind(-1.5, rep(1.5, 10)))
Domains <- rbind(Domains, cbind(-1.5, rep(1.5, 24 * 5)))
Domains <- rbind(Domains, cbind(0, rep(1, 24)))
fix_Delta_args <- fix_beta_args
fix_Delta_args$zeros <- c(1, 1)

Delta_list <- list(Delta = matrix(NA_real_, nrow = 5, ncol = 2), free = rep(TRUE, 10),
		fix_Delta_args = fix_Delta_args, 
		select = c(FALSE, rep(TRUE, 10), rep(FALSE, 144)))

res2_Harman74 <- new("restrictions.2ndorder", factors = c(5L, 2L),
			nvars = nrow(Domains), dof = 1L,
			Domains = Domains, model = "SEFA", method = "MLE",
			Xi = diag(c(0.5, 0.5)), Delta = Delta_list, Phi = diag(5),
			beta = beta_list, Theta2 = Theta2_list, 
			criteria = list(llik = FAiR:::FAiR_criterion_llik))

SEFA2_Harman74 <- Factanal(covmat = Harman74.cor, factors = 5, 
				restrictions = res2_Harman74)

## Compare CFA with example from library(sem)
semC <- rbind(
c(1.0000000, 0.8267444, 0.7745019, 0.4859589, 0.4635265, 0.4085056, 0.4732969, 0.4365315, 0.4264270),
c(0.8267444, 1.0000000, 0.7823011, 0.4908525, 0.4681942, 0.4126192, 0.4780630, 0.4409274, 0.4307211),
c(0.7745019, 0.7823011, 1.0000000, 0.4598352, 0.4386087, 0.3865455, 0.4478539, 0.4130649, 0.4035036),
c(0.4859589, 0.4908525, 0.4598352, 1.0000000, 0.6662540, 0.5871692, 0.4158054, 0.3835059, 0.3746288),
c(0.4635265, 0.4681942, 0.4386087, 0.6662540, 1.0000000, 0.5600648, 0.3966113, 0.3658028, 0.3573354),
c(0.4085056, 0.4126192, 0.3865455, 0.5871692, 0.5600648, 1.0000000, 0.3495332, 0.3223817, 0.3149195),
c(0.4732969, 0.4780630, 0.4478539, 0.4158054, 0.3966113, 0.3495332, 1.0000000, 0.5623102, 0.5492943),
c(0.4365315, 0.4409274, 0.4130649, 0.3835059, 0.3658028, 0.3223817, 0.5623102, 1.0000000, 0.5066255),
c(0.4264270, 0.4307211, 0.4035036, 0.3746288, 0.3573354, 0.3149195, 0.5492943, 0.5066255, 1.0000001))

R <- diag(9)
count <- 2
R[count, 1:(count-1)] <-    .828; count <- count + 1
R[count, 1:(count-1)] <-  c(.776, .779); count <- count + 1
R[count, 1:(count-1)] <-  c(.439, .493, .460); count <- count + 1
R[count, 1:(count-1)] <-  c(.432, .464, .425, .674); count <- count + 1
R[count, 1:(count-1)] <-  c(.447, .489, .443, .590, .541); count <- count + 1
R[count, 1:(count-1)] <-  c(.447, .432, .401, .381, .402, .288); count <- count + 1
R[count, 1:(count-1)] <-  c(.541, .537, .534, .350, .367, .320, .555); count <- count + 1
R[count, 1:(count-1)] <-  c(.380, .358, .359, .424, .446, .325, .598, .452)
R <- R + t(R)
diag(R) <- 1

varnames <- c('Sentences','Vocabulary','Sent.Completion','First.Letters','4.Letter.Words','Suffixes',
		'Letter.Series','Pedigrees', 'Letter.Group')
rownames(semC) <- colnames(semC) <- rownames(R) <- colnames(R) <- varnames
fixed <- matrix(0, nrow = 9, ncol = 3)
fixed[1:3,1] <- fixed[4:6,2] <- fixed[7:9,3] <- NA_real_

Domains <- rbind(cbind(-1, rep(1,3)), cbind(-1.5, rep(1.5, 9)), cbind(0, rep(1,9)))
Delta_list <- list(Delta = matrix(NA_real_, nrow = 3, ncol = 1), free = rep(TRUE, 3), num_free = 3,
			select = c(rep(TRUE, 3), rep(FALSE, nrow(Domains) - 3)))
beta_list <- list(beta = fixed, free = c(is.na(fixed)), num_free = sum(is.na(fixed)),
			select = c(rep(FALSE, 3), rep(TRUE, 9), rep(FALSE, 9)))
Theta2_list <- list(Theta2 = diag(9), select = c(rep(FALSE, 12), rep(TRUE, 9)))

res.Thur.1 <- new("restrictions.general", factors = c(3L, 1L),
			nvars = nrow(Domains), dof = 24L,
			Domains = Domains, model = "CFA", method = "MLE",
			Delta = Delta_list, Phi = diag(3),
			beta = beta_list, Theta2 = Theta2_list, 
			criteria = list(llik = FAiR:::FAiR_criterion_llik))

CFA <- Factanal(covmat = R, n.obs = 213, restrictions = res.Thur.1, factors = c(3,1), model = "CFA")
show(CFA)
summary(CFA)

fixed <- NA_real_ * fixed
beta_list <- list(beta = fixed, free = c(is.na(fixed)), num_free = sum(is.na(fixed)),
			select = c(rep(FALSE, 3), rep(TRUE, 27), rep(FALSE, 9)))
beta_list$fix_beta_args <- as.list(formals(FAiR:::FAiR_fix_coefficients))
beta_list$fix_beta_args$zeros <- rep(4,3)
Theta2_list <- list(Theta2 = diag(9), select = c(rep(FALSE, 30), rep(TRUE, 9)))
Domains <- rbind(cbind(-1, rep(1,3)), cbind(-1.5, rep(1.5, 27)), cbind(0, rep(1,9)))
res.Thur.2 <- new("restrictions.general", factors = c(3L, 1L),
			nvars = nrow(Domains), dof = 15L,
			Domains = Domains, model = "SEFA", method = "MLE",
			Delta = Delta_list, Phi = diag(3),
			beta = beta_list, Theta2 = Theta2_list, 
			criteria = list(llik = FAiR:::FAiR_criterion_llik))

SEFA <- Factanal(covmat = R, n.obs = 213, restrictions = res.Thur.2, factors = c(3,1), model = "SEFA")
show(SEFA)
summary(SEFA)
