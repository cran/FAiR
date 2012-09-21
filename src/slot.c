
/*  This file is part of FAiR, a program to conduct Factor Analysis in R
 *  Copyright (C) 2001-2005   The R Core Team.
 *  Copyright (C) 2012        Benjamin King Goodrich
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 3, or (at your option)
 *  any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  http://www.r-project.org/Licenses/
 */

#include <R.h>
#include <Rdefines.h>

SEXP FAiR_do_slot_assign(SEXP obj, SEXP name, SEXP value)
{
    return SET_SLOT(obj, name, value);
}



