/*
	Calculates the maximum value of the autocorrelation spectrum

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
 
#include "AC.h"

// ***
int AC(char *f, int FunctionsLength)
{
	int i=0, p=0, max=0;
	char *buf=NULL, *res=NULL;

	if(!f)
	{
		return -1;
	}

	buf=(char*)calloc(FunctionsLength,sizeof(char));
	res=(char*)calloc(FunctionsLength,sizeof(char));

	for(i=1;i<FunctionsLength;i++)
	{
		for(p=0;p<FunctionsLength;p++)
		{
			buf[p]=f[p^i];
		}

		res[i]=FunctionsLength-2*hd(buf,f,FunctionsLength);
	}

	for(max=0,i=0;i<FunctionsLength;i++)
	{
		if( max < abs(res[i]) )
		{
			max=abs(res[i]);
		}
	}

	if(buf)
		free(buf);

	if(res)
		free(res);

	return max;
}
// ***
