/*
	Calculates algebraic immunity of the vectorial boolean function given by 
	the component functions

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
 
#include "AlgebraicImmunity.h"

// ***
int algebraic_immunity(char **Functions, int deg, int FunctionsLength, int SboxBitIn, int SboxBitOut)
{
	int *tBuf=NULL, *cEq=NULL, x=0, y=0, z=0;
	int columns=1+(SboxBitIn+SboxBitOut), c=0, r=0, d=0, ret=0, counter=0, column=0;
	char **Matrix=NULL, **fs=NULL;

// columns in matrix
	for(c=2;c<=deg;c++)
	{
		columns+=(int)binom(SboxBitIn+SboxBitOut,c);
	}
// init functions
	tBuf=(int*)calloc((SboxBitIn+SboxBitOut)+2,sizeof(int)); // buffer for generating all combinations
	cEq=(int*)calloc((SboxBitIn+SboxBitOut),sizeof(int)); // current vector for combination
	fs=(char**)calloc(FunctionsLength,sizeof(char*)); // columns

	for(r=0;r<FunctionsLength;r++) // rows
	{
		fs[r]=(char*)calloc((SboxBitIn+SboxBitOut),sizeof(char));
	}

// filling SboxBitIn functions
	for(c=0;c<SboxBitIn;c++)
	{
		for(r=0;r<FunctionsLength;r++)
		{
			fs[r][c]=GET_BIT(r,SboxBitIn-1-c);
		}
	}

// filling SboxBitOut functions
	for(c=SboxBitIn;c<(SboxBitIn+SboxBitOut);c++)
	{
		for(r=0;r<FunctionsLength;r++)
		{
			fs[r][c]=Functions[c-SboxBitIn][r];
		}
	}

	init_matrix((char***)&Matrix,FunctionsLength,columns);

// filling matrix
	for(r=0;r<FunctionsLength;r++)
	{
		for(c=0;c<columns;c++)
		{
			Matrix[r][c]=1;
		}
	}

	for(d=1,column=1;d<=deg;d++)
	{
		inittwiddle(d, (SboxBitIn+SboxBitOut), tBuf);

		for( counter = 0; counter != (SboxBitIn+SboxBitOut)-d; counter++)
		{
			cEq[counter] = 0;
		}

		while( counter != (SboxBitIn+SboxBitOut) )
		{
			for(r=0;r<FunctionsLength;r++)
			{
				Matrix[r][column]&=fs[r][counter];
			}
			cEq[counter++] = 1;
		}

		column++;

		while(!twiddle(&x, &y, &z, tBuf))
		{
			cEq[x] = 1;
			cEq[y] = 0;

			for(counter = 0; counter != (SboxBitIn+SboxBitOut); counter++)
			{
				if(cEq[counter])
				{
					for(r=0;r<FunctionsLength;r++)
					{
						Matrix[r][column]&=fs[r][counter];
					}
				}
			} // counter
			column++;
		} // while
	} // d

// get rank matrix
	ret=rank_matrix(Matrix,FunctionsLength,columns);
	free_matrix(Matrix,FunctionsLength);
	
// delete temporary functions
	for(c=0;c<FunctionsLength;c++)
		if(fs[c])
			free(fs[c]);

	if(fs)
		free(fs);

	if(tBuf)
		free(tBuf);

	if(cEq)
		free(cEq);

	return (columns-ret) < 0 ? 0 : (columns-ret);
}
// ***
