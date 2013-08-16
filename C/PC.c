/*
	Calculates the value of propogation criterion

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
 
#include "PC.h"

// ***
int PC(int k,char *f, int FunctionsLength)
{
	unsigned int i=0, j=0;
	long long sum=0;
	char *g=NULL;
	int SboxBitOut=log2(FunctionsLength);

	if(!f)
	{
		return -1;
	}

	g=(char*)calloc(FunctionsLength,sizeof(char));

	for(i=0;i<FunctionsLength;i++)
	{
		if( (1 <= __builtin_popcountll(i)) && ( __builtin_popcountll(i) <= k ))
		{
			for(j=0;j<FunctionsLength;j++)
			{
				g[j]=f[j]^f[j^i];
			}

			sum=sumArray(g,FunctionsLength);

			if( sum != ((long long)1<<(SboxBitOut-1)) )
			{
				return 0;
			}
		}
	}

	if(g)
		free(g);

	return 1;
}
// ***