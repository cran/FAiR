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

## This file defines formal S4 classes

## NOTE: This file is meant to be read with 90 columns and 8 space tabs

## An object of class restrictions holds a lot of necessary information when using
## Factanal(). Several restrictions.* classes inherit from restrictions and cause
## Factanal() to estimate a somewhat different model.
setClass("restrictions", representation(factors   = "numeric",   # number of factors
					nvars     = "numeric",   # see genoud
					dof       = "integer",   # degrees of freedom
					Domains   = "matrix",    # bounds on parameters
					model     = "character", # what model?
					method    = "character"  # how estimated?
                                        ), package = "FAiR"
# 	showClass("restrictions")
# 	
# 	Slots:
# 									
# 	Name:    factors     nvars       dof   Domains     model    method
# 	Class:   numeric   numeric   integer    matrix character character
# 	
# 	Known Subclasses: "restrictions.factanal", "restrictions.orthonormal", 
# 	"restrictions.1storder", "restrictions.general", "restrictions.2ndorder"
)   ## end setClass("restrictions")

setClass("restrictions.factanal", representation("restrictions",
						    fast = "logical"), package = "FAiR"

)   ## end setClass("restrictions.mle")

setClass("restrictions.orthonormal", representation("restrictions",
							fast = "logical",
							Phi  = "matrix",
							beta = "list",
						      Theta2 = "list",
						    criteria = "list"), package = "FAiR"
# 	showClass("restrictions.factanal")
# 	
# 	Slots:
# 									
# 	Name:    factors     nvars       dof   Domains     model    method
# 	Class:   numeric   numeric   integer    matrix character character
# 	
# 	Extends: "restrictions"
)   ## end setClass("restrictions.orthonormal")

setClass("restrictions.1storder", representation("restrictions",
							Phi  = "matrix",
							beta = "list",
						      Theta2 = "list",
						    criteria = "list"), package = "FAiR"
# 	showClass("restrictions.1storder")
# 	
# 	Slots:
# 										
# 	Name:        Phi      beta    Theta2  criteria   factors     nvars       dof
# 	Class:    matrix      list      list      list   numeric   numeric   integer
# 					
# 	Name:    Domains     model    method
# 	Class:    matrix character character
# 	
# 	Extends: "restrictions"
)   ## end setClass("restrictions.1storder")

setClass("restrictions.general", representation("restrictions",
						      Delta  = "list",
							Phi  = "matrix",
							beta = "list",
						      Theta2 = "list",
						    criteria = "list"), package = "FAiR"
# 	showClass("restrictions.general")
# 	
# 	Slots:
# 										
# 	Name:      Delta       Phi      beta    Theta2  criteria   factors     nvars
# 	Class:      list    matrix      list      list      list   numeric   numeric
# 						
# 	Name:        dof   Domains     model    method
# 	Class:   integer    matrix character character
# 	
# 	Extends: "restrictions"
)   ## end setClass("restrictions.general")

setClass("restrictions.2ndorder", representation("restrictions",
						      Xi     = "matrix",
						      Delta  = "list",
							Phi  = "matrix",
							beta = "list",
						      Theta2 = "list",
						    criteria = "list"), package = "FAiR"
# 	showClass("restrictions.2ndorder")
# 	
# 	Slots:
# 										
# 	Name:         Xi     Delta       Phi      beta    Theta2  criteria   factors
# 	Class:    matrix      list    matrix      list      list      list   numeric
# 								
# 	Name:      nvars       dof   Domains     model    method
# 	Class:   numeric   integer    matrix character character
# 	
# 	Extends: "restrictions"
)   ## end setClass("restrictions.2ndorder")


## FA is the basic class for objects produced by Factanal() and Rotate(). It holds
## estimates and information about the completed estimation process. The two classes
## that extend FA deal with two-level factor analysis models.
setClass("FA", representation(  loadings     = "array",	    # variables x factors x 5
                                correlations = "array",     # factors   x factors x 3
                                trans_mats   = "array",     # factors   x factors x 2
                                uniquenesses = "numeric",   # vector of length(variables)
                                restrictions = "restrictions", # thing with many slots
                                vcov         = "matrix",    # parameters x parameters
                                zstats       = "list",      # length depends on model
                                scores       = "matrix",    # observations x factors
                                manifest     = "list",      # list containing LHS stuff
                                rgenoud      = "list",      # list produced by genoud()
                                model        = "character", # model estimated
                                method       = "character", # estimation method
                                call         = "language",  # call to function
                                seeds        = "matrix"     # PRNG seeds for genoud()
                             ), package = "FAiR"
# 	showClass("FA")
# 	
# 	Slots:
# 									
# 	Name:      loadings correlations   trans_mats uniquenesses restrictions
# 	Class:        array        array        array      numeric restrictions
# 									
# 	Name:          vcov       zstats       scores     manifest      rgenoud
# 	Class:       matrix         list       matrix         list         list
# 								
# 	Name:         model       method         call        seeds
# 	Class:    character    character     language       matrix
# 	
# 	Known Subclasses: "FA.general", "FA.2ndorder"
)   ## end setClass("FA")

setClass("FA.general", representation("FA",     restrictions = "restrictions.general",
						loadings_2nd = "matrix",
						uniquenesses_2nd = "numeric"), package = "FAiR"
# 	showClass("FA.general")
# 	
# 	Slots:
# 									
# 	Name:          restrictions         loadings_2nd     uniquenesses_2nd
# 	Class: restrictions.general               matrix              numeric
# 									
# 	Name:              loadings         correlations           trans_mats
# 	Class:                array                array                array
# 									
# 	Name:          uniquenesses                 vcov               zstats
# 	Class:              numeric               matrix                 list
# 									
# 	Name:                scores             manifest              rgenoud
# 	Class:               matrix                 list                 list
# 									
# 	Name:                 model               method                 call
# 	Class:            character            character             language
# 				
# 	Name:                 seeds
# 	Class:               matrix
# 	
# 	Extends: "FA"
)   ## end setClass("FA.general")

setClass("FA.2ndorder", representation("FA",    restrictions = "restrictions.2ndorder",
        					loadings_2nd = "array",
						uniquenesses_2nd = "numeric",
						correlations_2nd = "array"), package = "FAiR"
# 	showClass("FA.2ndorder")
# 	
# 	Slots:
# 										
# 	Name:           restrictions          loadings_2nd      uniquenesses_2nd
# 	Class: restrictions.2ndorder                 array               numeric
# 										
# 	Name:       correlations_2nd              loadings          correlations
# 	Class:                 array                 array                 array
# 										
# 	Name:             trans_mats          uniquenesses                  vcov
# 	Class:                 array               numeric                matrix
# 										
# 	Name:                 zstats                scores              manifest
# 	Class:                  list                matrix                  list
# 										
# 	Name:                rgenoud                 model                method
# 	Class:                  list             character             character
# 							
# 	Name:                   call                 seeds
# 	Class:              language                matrix
# 	
# 	Extends: "FA"
)   ## end setClass("FA.2ndorder")

## A class for the output produced by summary(FAobject)
setClass("summary.FA", representation( point_estimates = "list", 
						zstats = "list",
					  restrictions = "restrictions",
						call   = "call"), package = "FAiR"
# 	showClass("summary.FA")
# 	
# 	Slots:
# 									
# 	Name:  point_estimates          zstats    restrictions            call
# 	Class:            list            list    restrictions            call
)
