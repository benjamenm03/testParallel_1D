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
map<int, int> genArray(int iProc, int subSize, int size, int start, int end) {
    // pair<location on the x-axis, value>
    map<int, int> line;

    // fills the line with random numbers
    // each processor only fills their part of the line
    // rand() % (upper - lower + 1) + lower
    for (int i = start + (iProc * subSize); i < start + (iProc * subSize + subSize) && i < start + size; ++i) {
        line[i] = rand() % (end - start + 1) + start;
    }

    return line;
}

// fills an array with the processors rank for all the indices (correlating to the x-axis) 
// that the processor handles
vector<int> findLocations(map<int, int> *line, int size, int iProc) {
    // initialize an array with the size of the whole line filled with -1's
    vector<int> array (size, -1);

    // this array will be the result after combining all arrays from all processors
    vector<int> globalArray (size, -1);

    // fills in the array with iProc if the local processor is responsible for that region of the line
    // the rest of the line is still filled with -1's
    map<int, int>::iterator it = line->begin();
    while (it != line->end()) {
        array.insert(array.begin() + it->first, iProc);
        ++it;
    }

    // Combines all of the arrays that we just created and find the maximum values at every index across all arrays
    // and put put that max value into globalArray (this will get rid of the -1's and replace it with the correct 
    // processor number)
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Allreduce(array.data(), globalArray.data(), size, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

    return globalArray;
}



// Broadcast the location of one line (smaller line) (MPI_broadcast). 
// Then, create an array that is the same size. Array is filled with -1 on both processors. 
// Go through locations and ask, is this point on my processor? 
// If the point is on the processor, then fill in the rank of the processor onto the array. 
// The indices of the array are associated with the location of the line. 
// Do MPI_Allgather, take the max value, now every processor has a full array filled with 0, 1, 2, or 3, and there are no -1 left.
// Once you have the locations worked out, put a function on a processor that knows all of the data and try and get that data to all of the other processors

// All processors knows how big the line is (for example, it goes from x = 0 to x = 10)
// All processors knows that there are nProcs and nPts per processor
// Each processor then make their own line of nPts
// Each processor creates a line of nProcs * nPts and fill in the points that the local processor has
// Create a global grid of size nProcs * nPts that all processors know, and then try to share this grid

int main() {

    int nProcs, iProc;
    bool did_work;

    // NOTE - x cannot be negative

    // SET THE SIZE OF LINE 1
    // LINE WILL GO FROM x = 0 TO x = size - 1
    int start1 = 3;
    int end1 = 10;
    int size1 = end1 - start1 + 1;

    // SET THE SIZE OF LINE 2
    int start2 = 0;
    int end2 = 7;
    int size2 = end2 - start2 + 1;

    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &iProc);
    MPI_Comm_size(MPI_COMM_WORLD, &nProcs);

    // randomly calculate a seed for every processor
    int seed = rand() % (((iProc + 1) * 1000000) - (iProc * 1000000) + 1) + (iProc * 1000000);
    srand(seed);

    // calculates what processor takes care of how much of the array
    // to make more efficient, maybe split up the remainders to the processors so that it's more evenly distributed
    int subSize = size1 / nProcs;
    if (size1 % nProcs != 0) ++subSize;

    // prints out the processor rank and its corresponding line segment that it handles
    map<int, int> line1 = genArray(iProc, subSize, size1, start1, end1);
    for (int i = 0; i <= nProcs; ++i) {
        if (iProc == i) {
            cout << "Processor " << iProc << endl;
            map<int, int>::iterator it = line1.begin();
            while (it != line1.end()) {
                cout << "x = " << it->first << ", value = " << it->second << endl;
                ++it;
            }
            cout << endl;
        }
    }

    // finds what processor handles what part of the line and stores it in array
    // prints array
    vector<int> array1 = findLocations(&line1, size1, iProc);
    MPI_Barrier(MPI_COMM_WORLD);
    if (iProc == 0) {
        for (int i = start1; i < start1 + array1.size(); ++i) {
            cout << "x = " << i << ", processor = " << array1[i] << endl;
        }
    }


    MPI_Barrier(MPI_COMM_WORLD);
    // prints out the second line
    /*map<int, int> line2 = genArray(iProc, subSize, size2, start2, end2);
    if (iProc == 0) cout << "\n-------------------- LINE TWO --------------------\n" << endl;
    for (int i = 0; i <= nProcs; ++i) {
        if (iProc == i) {
            cout << "Processor " << iProc << endl;
            map<int, int>::iterator it = line2.begin();
            while (it != line2.end()) {
                cout << "x = " << it->first << ", value = " << it->second << endl;
                ++it;
            }
            cout << endl;
        }
    }*/

    // prints array for second line 
    /*vector<int> array2 = findLocations(&line2, size2, iProc);
    MPI_Barrier(MPI_COMM_WORLD);
    if (iProc == 0) {
        for (int i = 0; i < array2.size(); ++i) {
            cout << "x = " << i << ", processor = " << array2[i] << endl;
        }
    } */


    MPI_Finalize();
    return 0;
}