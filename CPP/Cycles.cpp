/*
	Finds cycles in a given array

AUTHORS:

- Maxim Storetvedt (2013-04-26): initial version

*/

/*****************************************************************************
 *       Copyright (C) 2013 Oleksandr Kazymyrov <oleksandr.kazymyrov@ii.uib.no>
*							Maxim Storetvedt <maxim.storetvedt@student.uib.no>
 *
 *  Distributed under the terms of the GNU General Public License (GPL)
 *  as published by the Free Software Foundation; either version 2 of
 *  the License, or (at your option) any later version.
 *                  http://www.gnu.org/licenses/
 ****************************************************************************/

#include "Cycles.h"

/*
 * Generates a random array and finds cycles in it using findCycles().
 */
// void testRun(){

// 	srand((unsigned)time(0));
// 	int n = (rand()%10);

// 	int nn = pow(2.0,n);
// 	//long long numbers[nn];

// 	//debug
// 	nn = 4;
// 	long long numbers[4] = {10, 3, 11, 6}; // {{0, 10}, {2, 11}, [{1, 3, 6}}
// 	//debug

// 	// for(int i = 0; i < nn; i++){
// 	//         numbers[i] = rand()%(nn-1);
// 	// }

// 	map<long long, vector<long long>*> cycles =     findCycles(numbers, (sizeof(numbers)/sizeof(numbers[0])));

// 	//Prints the cycle arrays.
// 	for(int i = 1; i < (int)cycles.size()+1; i++){
// 		if(!cycles[i]->empty())
// 			cout << "{";
// 			for(int j = 0; j < (int)cycles[i]->size(); j++){
// 				cout << cycles[i]->at(j);
// 				if(j < (int)cycles[i]->size()-1)
// 					cout << ",";
// 			}
// 			if(!cycles[i]->empty()){
// 				cout << "}";
// 				cout << endl;}
// 	}
// }

// int main() {

// 	testRun();

// 	return 0;
// }
/*
 *Method for joining together the corresponding values of two cycle arrays.
 */
bool joinCycle(vector<long long>* n, vector<long long>* f){

	int newEnd = n->at(n->size()-1); //End value of the new cycle array.
	int foundBeg = f->at(0); 	//Beginning value of the found cycle array.
	int foundEnd = f->at(f->size()-1); 	//End value of the found cycle array.

	//If the new array can naturally be added to the found one.
	if(newEnd == foundBeg){
		f->insert(f->begin(),n->begin(), n->end()-1);
		n->clear();

		return true;
	}

	//If the new array contains the sequences of the found one. (That is, the new
	//array contains the values, but not in the same order).
	else if (foundBeg == foundEnd){

		//Reorder and add new values.
		f->erase(f->begin());
		while(f->at(f->size()-1) != newEnd){
			f->push_back(f->at(0));
			f->erase(f->begin());
		}

		f->insert(f->begin(),n->begin(), n->end());
		n->clear();

		return true;
	}

	//If the new array only contains some the elements in the found array
	else{

		//Find the location that corresponds with the new array.
		int i = 0;
		int j = 0;
		bool foundI = false; //To stop adding to i after the location is found.
		while(f->at(j) != newEnd){

			if(f->at(i) == foundEnd)
				foundI = true;

			if(!foundI)
				i++;

			j++;
		}

		n->insert(n->end(),f->begin()+j+1, f->end());
		n->insert(n->end(),f->begin()+i+1, f->begin()+j+1);

		return false;
	}
}

/*
 *Find cycles in a given array.
 */
map<long long, vector<long long>*> findCycles(long long sbox[], long long s){

	map<long long, vector<long long>*> cycleArray; //We will store the found cycles here.
	vector<long long> cycleIndex; //keeps track of where each cycle is stored.

	cycleIndex.push_back(-1); //push one entry to the indexing array, as we will start at 1.

	int cycleNum = -1;
	int location = 0;
	int currentNum = 0;
	for(int i = 0; i < s; i++){

		if(sbox[i] >= 0){

			location = i;
			currentNum = 0;

			vector<long long>* tmp = new vector<long long>;
			cycleIndex.push_back(-1);

			//Iterate through the array
			while(true){

				//Update location and add value.
				currentNum = sbox[location];
				tmp->push_back(location);

				//If this position has previously not been visited.
				if (currentNum >= 0){

					//Set the position as the current cycle number.
					sbox[location] = cycleNum;

					//Check if the next location is outside the array. If it's not,
					//then continue iterating. If it is, add the value and stop iterating.
					if(currentNum < s)
						location = currentNum;
					else {
						tmp->push_back(currentNum);
						currentNum = cycleNum;
						break;}
				}

				//If this location is already visited.
				if (currentNum < 0){

					//Exit loop and add the found values.
					break;
				}
			}

			//If the location contains another cycle number than the current. (This indicates
			//that these numbers have already been added to a cycle array.
			if (currentNum != cycleNum){

				//Get the updated number.
				currentNum = cycleIndex[-(currentNum)];
				cycleIndex[-(cycleNum)] = currentNum;

				//Add the found values to the corresponding cycle. If not possible,
				//add the new cycle as its own entry in the cyclearray.
				if(!joinCycle(tmp,cycleArray[currentNum])){
					cycleArray[(long long) cycleArray.size() + 1] = tmp;
					cycleIndex[-(cycleNum)] = cycleArray.size();
				}

			}

			//If we've found an entirely new cycle, add it to the final cycleArray.
			else{

				cycleArray[cycleArray.size() + 1] = tmp;
				cycleIndex[-(cycleNum)] = cycleArray.size();
			}

			cycleNum--;
		}

	}

	return cycleArray;

}
