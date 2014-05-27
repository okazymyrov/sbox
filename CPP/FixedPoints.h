/*
	Finds fixed point in a given substitution

AUTHORS:

- Oleksandr Kazymyrov (2014-05-27): initial version

*/

/*****************************************************************************
 *       Copyright (C) 2013 Oleksandr Kazymyrov <oleksandr.kazymyrov@ii.uib.no>
 *
 *  Distributed under the terms of the GNU General Public License (GPL)
 *  as published by the Free Software Foundation; either version 2 of
 *  the License, or (at your option) any later version.
 *                  http://www.gnu.org/licenses/
 ****************************************************************************/

#include <cstdlib>
#include <iostream>
#include <vector>

using namespace std;

vector<unsigned long long> fixedPoints(unsigned long long[],unsigned long long);
