/*
	Auxiliary functions

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

#ifndef TOOLS_FUNCTIONS_H_
#define TOOLS_FUNCTIONS_H_

#include <stdlib.h>
#include <memory.h>
#include <stdio.h>
#include <math.h>

#define GET_BIT(x,b)	(((x)>>(b))&(1))
#define LOG()			printf("%d\n",__LINE__);

#ifdef __linux__
#define __rdtsc() \
       ({  unsigned int _a = 0 ; \
           unsigned int _d = 0; \
           asm volatile("rdtsc" : "=a" (_a), "=d" (_d)); \
           ((( unsigned long long)_a) | (((unsigned long long)_d) << 32)); })
#endif // __linux__

long long 			sumArray(char *arr, unsigned long long len);
unsigned long long		factorial(int n);
unsigned long long		binom(int m, int n);
int 				twiddle(int *x, int *y, int *z, int *p);
void 				inittwiddle(int m, int n, int *p);

#endif /* TOOLS_FUNCTIONS_H_ */
