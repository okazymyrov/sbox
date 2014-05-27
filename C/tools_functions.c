/*
	Auxiliary functions

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

#include "tools_functions.h"

long long sumArray(char *arr, unsigned long long len)
{
	unsigned long long i=0;
	unsigned long long sum=0;

	for(i=0;i<len;i++)
		sum+=arr[i];

	return sum;
}

int twiddle(int *x, int *y, int *z, int *p)
{
	register int i, j, k;
	j = 1;

	while(p[j] <= 0)
		j++;

	if(p[j-1] == 0)
	{
		for(i = j-1; i != 1; i--)
			p[i] = -1;

		p[j] = 0;
		*x = *z = 0;
		p[1] = 1;
		*y = j-1;
	}
	else
	{
		if(j > 1)
			p[j-1] = 0;
		do
		{
			j++;
		}
		while(p[j] > 0);

		k = j-1;
		i = j;

		while(p[i] == 0)
			p[i++] = -1;

		if(p[i] == -1)
		{
			p[i] = p[k];
			*z = p[k]-1;
			*x = i-1;
			*y = k-1;
			p[k] = -1;
		}
		else
		{
			if(i == p[0])
				return(1);
			else
			{
				p[j] = p[i];
				*z = p[i]-1;
				p[i] = 0;
				*x = j-1;
				*y = i-1;
			}
		}
	}

	return(0);
}

void inittwiddle(int m, int n, int *p)
{
	int i;
	p[0] = n+1;

	for(i = 1; i != n-m+1; i++)
		p[i] = 0;

	while(i != n+1)
	{
		p[i] = i+m-n;
		i++;
	}

	p[n+1] = -2;

	if(m == 0)
		p[1] = 1;
}

unsigned long long factorial(int n)
{
    int i;
    unsigned long long ret = 1;

    for (i = 1; i <= n; i++)
    {
        ret = ret * i;
    }

    return ret;
}

unsigned long long binom(int m, int n)
{
	return factorial(m)/(factorial(n)*factorial(m-n));
}
