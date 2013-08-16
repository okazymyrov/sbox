/*
Cycles one line description

AUTHORS:

- Oleksandr Kazymyrov (date in ISO year-month-day format): initial version

- person (date in ISO year-month-day format): short desc

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
#include <map>
#include <math.h>
using namespace std;

map<long long, vector<long long>*> findCycles(long long[],long long);
bool joinCycle (vector<long long>*,vector<long long>*);
