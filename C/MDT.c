/*
	Calculates the maximum value of differential table (\delta-uniformity)

AUTHORS:

- Oleksandr Kazymyrov (2013-08-03): initial version

*/

/*****************************************************************************
 *       Copyright (C) 2013 Oleksandr Kazymyrov <oleksandr.kazymyrov@ii.uib.no>
 *
 *  Distributed under the terms of the GNU General Public License (GPL)
 *  as published by the Free Software Foundation; either version 2 of
 *  the License, or (at your option) any later version.
 *                  http://www.gnu.org/licenses/
 ****************************************************************************/
 
#include "MDT.h"

// ***
unsigned long long d_uniform_xor_xor(unsigned long long *Sbox, unsigned long long SboxLength, unsigned long long TableSize, unsigned long long delta)
{
	unsigned long long max=0, j=0, i=0;
	unsigned long long *DifTable=NULL;

	DifTable=(unsigned long long*)calloc(TableSize,sizeof(unsigned long long));

	for(j=1;j<SboxLength;j++)
	{
		for(i=0;i<SboxLength;i++)
			DifTable[Sbox[i^j]^Sbox[i]]++;

		for(i=0;i<TableSize;i++)
		{
			if(DifTable[i]>delta)
			{
				if(DifTable)
					free(DifTable);
				return 0;
			}

			if(max<DifTable[i])
				max=DifTable[i];

			DifTable[i]=0;
		}
	}

	if(DifTable)
		free(DifTable);
	
	return max;
}
// ***
unsigned long long d_uniform_add_xor(unsigned long long *Sbox, unsigned long long SboxLength, unsigned long long TableSize, unsigned long long n, unsigned long long m, unsigned long long delta)
{
	unsigned long long max=0, j=0, i=0, modulo=((unsigned long long)1<<n)-1;
	unsigned long long *DifTable = NULL;

	DifTable = (unsigned long long*)calloc(TableSize,sizeof(unsigned long long));

	for(j=1;j<SboxLength;j++)
	{
		for(i=0;i<SboxLength;i++)
			DifTable[Sbox[(i+j)&modulo]^Sbox[i]]++;

		for(i=0;i<TableSize;i++)
		{
			if(DifTable[i]>delta)
			{
				if(DifTable)
					free(DifTable);
				return 0;
			}

			if(max<DifTable[i])
				max=DifTable[i];

			DifTable[i]=0;
		}
	}

	if(DifTable)
		free(DifTable);
	
	return max;
}
// ***
unsigned long long d_uniform_xor_add(unsigned long long *Sbox, unsigned long long SboxLength, unsigned long long TableSize, unsigned long long n, unsigned long long m, unsigned long long delta)
{
	unsigned long long max=0, j=0, i=0, modulo=((unsigned long long)1<<m)-1;
	unsigned long long *DifTable = NULL;

	DifTable = (unsigned long long*)calloc(TableSize,sizeof(unsigned long long));

	for(j=1;j<SboxLength;j++)
	{
		for(i=0;i<SboxLength;i++)
			DifTable[(Sbox[i^j] - Sbox[i])&modulo]++;

		for(i=0;i<TableSize;i++)
		{
			if(DifTable[i]>delta)
			{
				if(DifTable)
					free(DifTable);
				return 0;
			}

			if(max<DifTable[i])
				max=DifTable[i];

			DifTable[i]=0;
		}
	}

	if(DifTable)
		free(DifTable);
	
	return max;
}
// ***
unsigned long long d_uniform_add_add(unsigned long long *Sbox, unsigned long long SboxLength, unsigned long long TableSize, unsigned long long n, unsigned long long m, unsigned long long delta)
{
	unsigned long long max=0, j=0, i=0, modulo_in=((unsigned long long)1<<n)-1, modulo_out=((unsigned long long)1<<m)-1;
	unsigned long long *DifTable = NULL;

	DifTable = (unsigned long long*)calloc(TableSize,sizeof(unsigned long long));

	for(j=1;j<SboxLength;j++)
	{
		for(i=0;i<SboxLength;i++)
			DifTable[(Sbox[(i+j)&modulo_in] - Sbox[i])&modulo_out]++;

		for(i=0;i<TableSize;i++)
		{
			if(DifTable[i]>delta)
			{
				if(DifTable)
					free(DifTable);
				return 0;
			}

			if(max<DifTable[i])
				max=DifTable[i];

			DifTable[i]=0;
		}
	}

	if(DifTable)
		free(DifTable);
	
	return max;
}
// ***