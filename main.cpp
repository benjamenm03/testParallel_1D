// main.cpp file for testing one-dimensional examples of MPI programs for Aether

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <map>
#include <unistd.h>

using namespace std;


// creates an array filled with random numbers
map<int, int> genArray(int size, int range_start, int range_end, int seed) {
    // pair<location on the x-axis, value>
    map<int, int> line;
    int nProcs, iProc;
    MPI_Comm_rank(MPI_COMM_WORLD, &iProc);
    MPI_Comm_size(MPI_COMM_WORLD, &nProcs);

    // calculates what processor takes care of how much of the array
    // to make more efficient, maybe split up the remainders to the processors so that it's more evenly distributed
    int subSize = size / nProcs;
    if (size % nProcs != 0) ++subSize;

    // set the seed for the random number generator
    srand(seed);

    // fills the line with random numbers
    // each processor only fills their part of the line
    // rand() % (upper - lower + 1) + lower
    for (int i = iProc * subSize; i < iProc * subSize + subSize && i < size; ++i) {
        line[i] = rand() % (range_end - range_start + 1) + range_start;
    }

    return line;
}

// fills an array with the processors rank for all the indices (correlating to the x-axis) 
// that the processor handles
void findLocations(vector<int> *array, map<int, int> *line) {
    int iProc;
    MPI_Comm_rank(MPI_COMM_WORLD, &iProc);

    map<int, int>::iterator it = line->begin();
    while (it != line->end()) {
        array->insert(array->begin() + it->first, iProc);
        ++it;
    }
}

// Broadcast the location of one line (smaller line) (MPI_broadcast). 
// Then, create an array that is the same size. Array is filled with -1 on both processors. 
// Go through locations and ask, is this point on my processor? 
// If the point is on the processor, then fill in the rank of the processor onto the array. 
// The indices of the array are associated with the location of the line. 
// Do MPI_Allgather, take the max value, now every processor has a full array filled with 0, 1, 2, or 3, and there are no -1 left.
// // OR have a global array that can pass by reference and have all the processors fill in with????
// Once you have the locations worked out, put a function on a processor that knows all of the data and try and get that data to all of the other processors


// I NEED TO FIGURE OUT HOW TO STORE THE LINE. RIGHT NOW, THE MAP GETS DESTROYED AFTER genArray FINISHES
int main() {

    int nProcs, iProc;
    bool did_work;

    // SET THE SIZE OF THE LINE
    // LINE WILL GO FROM x = 0 TO x = size - 1
    int size = 10;
    vector<int> locations;
    locations.resize(size);

    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &iProc);
    MPI_Comm_size(MPI_COMM_WORLD, &nProcs);

    // prints out the processor rank and its corresponding line segment that it handles
    map<int, int> line = genArray(size, 0, 100, 123);
    for (int i = 0; i <= nProcs; ++i) {
        if (iProc == i) {
            cout << "Processor " << iProc << endl;
            map<int, int>::iterator it = line.begin();
            while (it != line.end()) {
                cout << "x = " << it->first << ", value = " << it->second << endl;
                ++it;
            }
            cout << endl;
        }
    }

    findLocations(&locations, &line);

    MPI_Finalize();
    return 0;
}