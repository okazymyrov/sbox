/*
	Auxiliary functions to work with matrices

AUTHORS:

- Oleksandr Kazymyrov (2013-08-03): initial version

*/

/*****************************************************************************
 *       Copyright (C) 2013 Oleksandr Kazymyrov <okazymyrov@gmail.com>
 *
 *  Distributed under the terms of the GNU General Public License (GPL)
 *  as published by the Free Software Foundation; either version 2 of
 *  the License, or (at your option) any later version.
 *                  http://www.gnu.org/licenses/
 ****************************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "tools_functions.h"

#define min(a,b)    (((a) < (b)) ? (a) : (b))

int					init_matrix(char ***m, int rows, int columns);
int					free_matrix(char **m, int rows);

int					xor_vectos(char *a,char *b,char *d,int length);
int					swap_matrix_row(char **m,int r1, int r2);
int					row_reduce_matrix(char **m, int rows, int columns);
int					rank_matrix(char **m, int rows, int columns);
