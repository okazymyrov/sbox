/*
    Calculates the maximum value of linear approximation table

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
 
#include "MLT.h"

// Function's names have gotten from article "Fast Computation of Approximation Tables" (Krzysztof Chmiel)
inline unsigned long long BIT_XOR(unsigned long long value, unsigned long long mask)
{
    value = value & mask;
    unsigned long long xorbit = value;
    unsigned long long i=sizeof(unsigned long long)<<2;

    while(i != 0)
    {
	xorbit = (xorbit&(((unsigned long long)1<<i)-1)) ^ (xorbit>>i);
	i>>=1;
    }
    
    return xorbit;
}

inline void TP(long long *v0, long long *v1)
{
    if(*v0 == 0 )
    {
	if(*v1 == 0)
	{
    	    *v0 = 1;
    	    *v1 = 0;
    	}
    }
    else
    {
	if(*v1 == 0)
	{
    	    *v0 = 0;
    	    *v1 = -1;
    	}
	else
	{
    	    *v0 = -1;
    	    *v1 = 0;
    	}
    }
}

inline void INI_TA(unsigned long long *sbox, long long InLength, long long OutLength, long long* col, unsigned long long mask_y)
{
    unsigned long long mask_x;
    for(mask_x = 0; mask_x < InLength; mask_x++)
    {
        col[mask_x] = BIT_XOR(sbox[mask_x], mask_y); 
    }

    for(mask_x = 0; mask_x < InLength - 1; mask_x += 2)
    {
        TP(&col[mask_x], &col[mask_x + 1]);
    }
}

inline void SUMSUB_TA(long long *col, long long i1, long long j1, long long i2, long long j2)
{
    static long long temp; 
    long long i=0;
    for(i = 0; i <= (j1 - i1); i++)
    {
        temp = col[i1 + i];

        col[i1 + i] = temp + col[i2 + i];
        col[i2 + i] = temp - col[i2 + i];
    }
}

void CALC_TA(long long *col, size_t len, long long i, long long j)
{
    long long k;
    if(j - i > 2)
    {
        k = (i + j) >> 1;
        CALC_TA(col, len, i, k);
        CALC_TA(col, len, k+1, j);
        SUMSUB_TA(col, i, k, k+1, j);
    }
}

unsigned long long c_lambda(unsigned long long *sbox, long long InBit, long long OutBit)
{
    long long *col = NULL;
    unsigned long long max=0, i=0, mask_y=0;
    long long InLength=1<<InBit, OutLength=1<<OutBit;

    col=(long long *)calloc(InLength,sizeof(long long));

    INI_TA(sbox, InLength, OutLength, col, 0);
    CALC_TA(col, InLength, 0, InLength-1);

    for(i = 1; i < InLength; i++)
    {
        if (abs(col[i]) > max)
            max=abs(col[i]);
    }

    for(mask_y = 1; mask_y < OutLength; mask_y++)
    {
        INI_TA(sbox, InLength, OutLength, col, mask_y);
        CALC_TA(col, InLength, 0, InLength-1);

        for(i = 0; i < InLength; i++)
        {
            if ( abs(col[i]) > max)
                max=abs(col[i]);
        }
    }

    if(col)
        free(col);

    return max;
}
// End lin_Sbox
