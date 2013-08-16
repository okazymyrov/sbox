/*
	Calculates algebraic immunity of the vectorial boolean function given by 
	the component functions

AUTHORS:

- Oleksandr Kazymyrov (date in ISO year-month-day format): initial version

*/

/*****************************************************************************
 *       Copyright (C) 2013 Oleksandr Kazymyrov <oleksandr.kazymyrov@ii.uib.no>
 *
 *  Distributed under the terms of the GNU General Public License (GPL)
 *  as published by the Free Software Foundation; either version 2 of
 *  the License, or (at your option) any later version.
 *                  http://www.gnu.org/licenses/
 ****************************************************************************/

#include <stdio.h>

#include "tools_functions.h"

int algebraic_immunity(char **Functions, int deg, int FunctionsLength, int SboxBitIn, int SboxBitOut);