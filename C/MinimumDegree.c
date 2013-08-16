/*
	Calculates the minimum degree of the substitution

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
 
#include "MinimumDegree.h"

// ***
int ANF(char *f, char **anf, int FunctionsLength)
{
	int i=0, j=0;
	char *t_f=NULL, *u=NULL, *v=NULL;
	int SboxBitOut=log2(FunctionsLength);

	if(!f || !anf)
	{
		return -1;
	}

	t_f=(char *)calloc(FunctionsLength,sizeof(char));
	u=(char *)calloc(FunctionsLength>>1,sizeof(char));
	v=(char *)calloc(FunctionsLength>>1,sizeof(char));

	memcpy(t_f,f,FunctionsLength);

	for( i = 0; i < SboxBitOut; i++)
	{
		for( j = 0; j < (FunctionsLength>>1); j++)
		{
			v[j] = t_f[ 2*j ];
			u[j] = t_f[ 2*j ] ^ t_f[ 2*j + 1 ];
		}

		for( j = 0; j < (FunctionsLength>>1); j++)
		{
			t_f[j] = v[j];
			t_f[j+(FunctionsLength>>1)] = u[j];
		}
	}

	memcpy(*anf,t_f,FunctionsLength);

	if(u)
		free(u);
	if(v)
		free(v);
	if(t_f)
		free(t_f);

	return 1;
}
// ***
int minimum_degree(char *f, int FunctionsLength)
{
	int i=0, f_deg=0;
	char *anf=NULL;

	anf=(char*)calloc(FunctionsLength,sizeof(char));

	ANF(f,&anf,FunctionsLength);

	for(f_deg=0,i=FunctionsLength-1;i>=0;i--)
	{
		if( anf[i] )
			if( f_deg<__builtin_popcountll(i) )
				f_deg=__builtin_popcountll(i);
	}

	if(anf)
		free(anf);

	return f_deg;
}
// ***