/*
	Calculates tthe correlation immunity of the boolean function

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
 
#include "CI.h"

// ***
int hd(char *f, char *g, int FunctionsLength)
{
	int ret=0, i=0;

	for(ret=0,i=0;i<FunctionsLength;i++)
		if(f[i]^g[i])
			ret++;

	return ret;
}
// ***
int WHT(int w, char *f, int FunctionsLength)
{
	int i=0, p=0, ret=0;
	char *buf=NULL;
	int SboxBitOut=log2(FunctionsLength);

	buf=(char*)calloc(FunctionsLength,sizeof(char));

	for(p=0;p<FunctionsLength;p++)
	{
		buf[p]=0;

		for(i=0;i<SboxBitOut;i++)
		{
			buf[p]^=(GET_BIT(w,i)&GET_BIT(p,i));
		}		
	}

	ret=FunctionsLength-2*hd(buf,f,0);

	if(buf)
		free(buf);

	return ret;
}
// ***
int CI(int m,char *f, int FunctionsLength)
{
	int w=0;

	if(!f)
	{
		return -1;
	}

	for(w=0;w<FunctionsLength;w++)
	{
		if( (1 <= __builtin_popcountll(w)) || (__builtin_popcountll(w) <= m) )
		{
			if( WHT(w,f,FunctionsLength) != 0 )
			{
				return 0;
			}
		}
	}

	return 1;
}
// ***