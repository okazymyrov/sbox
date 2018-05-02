/*
	Finds cycles in a given array

AUTHORS:

- Maxim Storetvedt (2013-04-26): initial version

*/

/*****************************************************************************
 *       Copyright (C) 2013 Oleksandr Kazymyrov <okazymyrov@gmail.com>
 *			    Maxim Storetvedt <maxim.storetvedt@student.uib.no>
 *
 *  Distributed under the terms of the GNU General Public License (GPL)
 *  as published by the Free Software Foundation; either version 2 of
 *  the License, or (at your option) any later version.
 *                  http://www.gnu.org/licenses/
 ****************************************************************************/

#include <cstdlib>
#include <iostream>
#include <vector>
#include <map>
#include <math.h>
using namespace std;

map<long long, vector<long long>*> findCycles(long long[],long long);
bool joinCycle (vector<long long>*,vector<long long>*);
