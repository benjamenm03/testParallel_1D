// main.cpp file for testing one-dimensional examples of MPI programs for Aether

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <map>
#include <unistd.h>

using namespace std;

/*void mapSum(void *inputBuffer, void *outputBuffer, int *length, MPI_Datatype *datatype) {
    map<int, int> *in = static_cast<map<int, int> *> (inputBuffer);
    map<int, int> *out = static_cast<map<int, int> *> (outputBuffer);
    
    map<int, int>::iterator input = in->begin();
    map<int, int>::iterator output = out->begin();

    // CANT ACCESS THE MAPS FOR SOME REASON
    while (input != in->end()) {
        cout << "HELP x = " << input->first << ", value = " << input->second << endl;
        ++input;
    }

    for (int i = 0; i < *length; ++i) {
        cout << "HELP x = " << input->first << ", value = " << input->second << endl;
        if (i < *length - 1) ++input;
    }

    return;
    // SOMETHING WRONG WITH THIS
    while (input != in->end() && output != out->end()) {
        if (input->second > output->second) output->second = input->second;
        ++input;
        ++output;
    }
}*/

// creates an array filled with random numbers
map<int, int> genArray(int iProc, int subSize, int size, int start, int end) {
    // pair<location on the x-axis, value>
    map<int, int> grid;

    // fills the grid with random numbers
    // each processor only fills their part of the grid
    // rand() % (upper - lower + 1) + lower
    for (int i = start + (iProc * subSize); i < start + (iProc * subSize + subSize) && i < start + size; ++i) {
        grid[i] = rand() % (end - start + 1) + start;
    }

    return grid;
}

// fills an array with the processors rank for all the indices (correlating to the x-axis) 
// that the processor handles
map<int, int> findLocations(map<int, int> *grid, int size, int start, int iProc) {
    // initialize a vector with the size of the whole grid filled with -1's
    vector<int> array (size, -1);

    // this vector will be the result after combining all arrays from all processors
    vector<int> globalArray (size, -1);

    // fills in the vector with iProc if the local processor is responsible for that region of the grid
    // the rest of the grid is still filled with -1's
    map<int, int>::iterator it = grid->begin();
    while (it != grid->end()) {
        array[it->first - start] = iProc;
        ++it;
    }

    // Combines all of the vectors that we just created and find the maximum values at every index across all vectors
    // and put that max value into globalArray (this will get rid of the -1's and replace it with the correct 
    // processor rank)
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Allreduce(array.data(), globalArray.data(), array.size(), MPI_INT, MPI_MAX, MPI_COMM_WORLD);

    // turns the vector back into a map so that we have the right x values
    map<int, int> answer;
    for (int i = 0; i < globalArray.size(); ++i) {
        answer[i + start] = globalArray[i];
    }

    return answer;
}



// Broadcast the location of one grid (smaller grid) (MPI_broadcast). 
// Then, create an array that is the same size. Array is filled with -1 on both processors. 
// Go through locations and ask, is this point on my processor? 
// If the point is on the processor, then fill in the rank of the processor onto the array. 
// The indices of the array are associated with the location of the grid. 
// Do MPI_Allgather, take the max value, now every processor has a full array filled with 0, 1, 2, or 3, and there are no -1 left.
// Once you have the locations worked out, put a function on a processor that knows all of the data and try and get that data to all of the other processors

// All processors knows how big the grid is (for example, it goes from x = 0 to x = 10)
// All processors knows that there are nProcs and nPts per processor
// Each processor then make their own grid of nPts
// Each processor creates a grid of nProcs * nPts and fill in the points that the local processor has
// Create a global grid of size nProcs * nPts that all processors know, and then try to share this grid

int main() {

    int nProcs, iProc;
    bool did_work;

    // SET THE SIZE OF grid 1
    int start1 = -2;
    int end1 = 5;
    int size1 = end1 - start1 + 1;

    // SET THE SIZE OF grid 2
    int start2 = 1;
    int end2 = 8;
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

    // prints out the processor rank and its corresponding grid segment that it handles
    map<int, int> grid1 = genArray(iProc, subSize, size1, start1, end1);
    for (int i = 0; i <= nProcs; ++i) {
        if (iProc == i) {
            cout << "Processor " << iProc << endl;
            map<int, int>::iterator it = grid1.begin();
            while (it != grid1.end()) {
                cout << "x = " << it->first << ", value = " << it->second << endl;
                ++it;
            }
            cout << endl;
        }
    }

    // finds what processor handles what part of the grid and stores it in array
    // prints array
    MPI_Barrier(MPI_COMM_WORLD);
    map<int, int> array1 = findLocations(&grid1, size1, start1, iProc);
    MPI_Barrier(MPI_COMM_WORLD);
    if (iProc == 0) {
        map<int, int>::iterator it = array1.begin();
        while (it != array1.end()) {
            cout << "x = " << it->first << ", value = " << it->second << endl;
            ++it;
        }
    }


    MPI_Barrier(MPI_COMM_WORLD);
    // prints out the second grid
    map<int, int> grid2 = genArray(iProc, subSize, size2, start2, end2);
    if (iProc == 0) cout << "\n-------------------- grid TWO --------------------\n" << endl;
    for (int i = 0; i <= nProcs; ++i) {
        if (iProc == i) {
            cout << "Processor " << iProc << endl;
            map<int, int>::iterator it = grid2.begin();
            while (it != grid2.end()) {
                cout << "x = " << it->first << ", value = " << it->second << endl;
                ++it;
            }
            cout << endl;
        }
    }

    // prints array for second grid 
    MPI_Barrier(MPI_COMM_WORLD);
    map<int, int> array2 = findLocations(&grid2, size2, start2, iProc);
    MPI_Barrier(MPI_COMM_WORLD);
    if (iProc == 0) {
        map<int, int>::iterator it = array2.begin();
        while (it != array2.end()) {
            cout << "x = " << it->first << ", value = " << it->second << endl;
            ++it;
        }
    }


    MPI_Finalize();
    return 0;
}