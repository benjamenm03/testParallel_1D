#include <vector>
#include <map>
#include <iostream>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

using namespace std;

/// @brief creates a map that goes from x = start to x = end (inclusive) 
///        with each x-value getting a randomly generated double
/// @param iProc the rank of the processor that's calling the function
/// @param subSize the size of the grid that each processor is responsible for
/// @param size the total size of the grid
/// @param start the starting x-value of the grid
/// @param end the ending x-value of the grid
/// @return returns the map
map<double, double> genArray(int iProc, int subSize, int size, int start, int end, int begin, int finish) {
    // pair<location on the x-axis, value>
    map<double, double> grid;

    // fills the grid with random numbers
    // each processor only fills their part of the grid
    // rand() % (upper - lower + 1) + lower
    for (int i = start + (iProc * subSize); i < start + (iProc * subSize + subSize) && i < start + size; ++i) {
        grid[i] = random() % (begin - finish + 1) + begin;
    }

    return grid;
}



/// @brief fills an array with the processors rank for all the indices 
///        (correlating to the x-axis) that the processor handles
/// @param grid the grid that each processor has
/// @param size the total size of the grid
/// @param start the starting x-value of the grid
/// @param iProc the processor's rank that's calling the function
/// @return a map that contains information about what processor knows what
map<double, double> findLocations(map<double, double> *grid, int size, int start, int iProc) {
    // initialize a vector with the size of the whole grid filled with -1's
    vector<double> array (size, -1);

    // this vector will be the result after combining all arrays from all processors
    vector<double> globalArray (size, -1);

    // fills in the vector with iProc if the local processor is responsible for that region of the grid
    // the rest of the grid is still filled with -1's
    map<double, double>::iterator it = grid->begin();
    while (it != grid->end()) {
        array[it->first - start] = iProc;
        ++it;
    }

    // Combines all of the vectors that we just created and find the maximum values at every index across all vectors
    // and put that max value into globalArray (this will get rid of the -1's and replace it with the correct 
    // processor rank)
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Allreduce(array.data(), globalArray.data(), array.size(), MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    // turns the vector back into a map so that we have the right x values
    map<double, double> answer;
    for (int i = 0; i < globalArray.size(); ++i) {
        answer[i + start] = globalArray[i];
    }

    return answer;
}



/// @brief prints the whole grid
/// @param array the map that contains information about what processor knows what
/// @param grid the grid that each processor has
/// @param iProc the processor rank that's calling the function
void printGrid(map<double, double> *array, map<double, double> *grid, int iProc) {
    map<double, double>::iterator it = array->begin();

    while (it != array->end()) {
        if (iProc == it->second) {
            cout << "x = " << it->first << ", value = " << grid->at(it->first) << endl;
        }
        ++it;
    }
}



/// @brief find the rank of the processor that contains the data at x = xPos
/// @param array the array that contains information about what processor knows what
/// @param xPos the x position that we are interested in
/// @return the rank of the processor tha contains the data at x = xPos.
///         If the xPos cannot be found, return -1
int findProc(map<double, double> *array, double xPos) {
    map<double, double>::iterator it = array->begin();
    
    while (it != array->end()) {
        if (it->first == xPos) return it->second;
        ++it;
    }

    return -1;
}



/// @brief retrieves data from another grid (MPI message passing)
/// @param receive the map that contains information on what processor has what data for the grid that wants the data
/// @param send the map that contains information on what processor has what data for the grid that needs to send the data
/// @param grid the grid that needs to send the data
/// @param iProc the rank of the processor that's executing the function
/// @param xPos the x-position that we are interested in getting data to and from
/// @param answer the data will be stored in this variable on the processor that wanted the data
/// @return 
int getValue(map<double, double> *receive, map<double, double> *send, map<double, double> *grid, 
              int iProc, double xPos, double *answer) {
    int rcv = findProc(receive, xPos);
    int snd = findProc(send, xPos);

    if (iProc == snd) {
        MPI_Send(&grid->at(xPos), 1, MPI_DOUBLE, rcv, 0, MPI_COMM_WORLD);
    }
    else if (iProc == rcv) {
        MPI_Recv(answer, 1, MPI_DOUBLE, snd, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    return rcv;
}




map<double, double> getCoeff(map<double, double> *receive, map<double, double> *get) {

}

