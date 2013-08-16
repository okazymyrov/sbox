/*
	Calculates the sum-of-squares indicator

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

#include "SSI.h"

// ***
unsigned long long SumOfSquares(char *f, int FunctionsLength)
{
	unsigned long long ret=0;
	char *buf=NULL, *res=NULL;
	int i=0, p=0;

	buf=(char*)calloc(FunctionsLength,sizeof(char));
	res=(char*)calloc(FunctionsLength,sizeof(char));

	for(i=1;i<FunctionsLength;i++)
	{
		for(p=0;p<FunctionsLength;p++)
		{
			buf[p]=1-2*(f[p]^f[p^i]);
		}

		res[i]=sumArray(buf,FunctionsLength);
	}

	for(ret=0,i=0;i<FunctionsLength;i++)
	{
		ret+=res[i]*res[i];
	}

	ret+=FunctionsLength*FunctionsLength;

	if(buf)
		free(buf);

	if(res)
		free(res);

	return ret;
}
// ***