/*
	Calculates the maximum value of differential table (\delta-uniformity)

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

unsigned long long d_uniform(unsigned long long *Sbox, unsigned long long SboxLength, unsigned long long TableSize, unsigned long long delta);