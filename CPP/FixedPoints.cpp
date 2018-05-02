/*
	Finds fixed point in a given substitution

AUTHORS:

- Oleksandr Kazymyrov (2014-05-27): initial version

*/

/*****************************************************************************
 *       Copyright (C) 2013 Oleksandr Kazymyrov <okazymyrov@gmail.com>
 *
 *  Distributed under the terms of the GNU General Public License (GPL)
 *  as published by the Free Software Foundation; either version 2 of
 *  the License, or (at your option) any later version.
 *                  http://www.gnu.org/licenses/
 ****************************************************************************/

#include "FixedPoints.h"

/*
 *Find cycles in a given array.
 */
vector<unsigned long long> fixedPoints(unsigned long long sbox[], unsigned long long length)
{
	vector<unsigned long long> fPoints;

	for(unsigned long long i = 0; i < length; i++)
	{
		if(sbox[i] == i)
		{
			fPoints.push_back(i);
		}
	}

	return fPoints;
}
