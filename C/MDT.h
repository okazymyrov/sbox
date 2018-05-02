/*
	Calculates the maximum value of differential table (\delta-uniformity)

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

unsigned long long d_uniform_xor_xor(unsigned long long *Sbox, unsigned long long SboxLength, unsigned long long TableSize, unsigned long long delta);
unsigned long long d_uniform_add_xor(unsigned long long *Sbox, unsigned long long SboxLength, unsigned long long TableSize, unsigned long long n, unsigned long long m, unsigned long long delta);
unsigned long long d_uniform_xor_add(unsigned long long *Sbox, unsigned long long SboxLength, unsigned long long TableSize, unsigned long long n, unsigned long long m, unsigned long long delta);
unsigned long long d_uniform_add_add(unsigned long long *Sbox, unsigned long long SboxLength, unsigned long long TableSize, unsigned long long n, unsigned long long m, unsigned long long delta);