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
map<double, double> genArray(int iProc, int subSize, int size, double interval, int start, int end, int begin, int finish) {
    // pair<location on the x-axis, value>
    map<double, double> grid;

    // fills the grid with random numbers
    // each processor only fills their part of the grid
    // rand() % (upper - lower + 1) + lower
    for (double i = start + (iProc * subSize); i < start + (iProc * subSize + subSize) && i < start + size; i += interval) {
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
map<double, double> findLocations(map<double, double> *grid, int size, double interval, int start, int iProc) {
    // initialize a vector with the size of the whole grid filled with -1's
    vector<double> array (size / interval, -1);

    // this vector will be the result after combining all arrays from all processors
    vector<double> globalArray (size / interval, -1);

    // fills in the vector with iProc if the local processor is responsible for that region of the grid
    // the rest of the grid is still filled with -1's

    // THIS PART IS WRONG - VECTORS CANNOT HAVE INTERVALS OTHER THAN 1
    map<double, double>::iterator it = grid->begin();
    while (it != grid->end()) {
        array[abs(it->first - start) / interval] = iProc;
        ++it;
    }
    

    // Combines all of the vectors that we just created and find the maximum values at every index across all vectors
    // and put that max value into globalArray (this will get rid of the -1's and replace it with the correct 
    // processor rank)
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Allreduce(array.data(), globalArray.data(), array.size(), MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    // turns the vector back into a map so that we have the right x values
    map<double, double> answer;
    vector<double>::iterator it2 = globalArray.begin();
    for (double i = start; i < start + size; i += interval) {
        answer[i] = *it2;
        ++it2;
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




// The coefficient is always with respect to the point that is smaller than the point we're interested in
// ex: x = 2, x = 3. If we want to find x = 2.4, then coeff = 40

/// @brief find the coefficient for the receive grid with respect to the get grid
/// @param receive the ownership map for the grid we want to get coefficients for (that way we have access to all indices)
/// @param get the ownership map that we are using to calculate the coefficients
/// @param interval the interval of the receive grid (spacing between each x-value)
/// @return returns the coefficient map
map<double, double> findCoeff(map<double, double> *receive, map<double, double> *get) {
    int size = receive->size();
    map<double, double> coeff;
    map<double, double>::iterator recv = receive->begin();

    // iterate through the receive array (the coefficient array is going to contain the same x-values as the receive array)
    while (recv != receive->end()) {

        // if a x-value goes above or below the bounds of the get array, set the coefficient as -1
        if (recv->first < get->begin()->first || recv->first > get->rbegin()->first) {
            coeff[recv->first] = -1;
            ++recv;
            continue;
        }

        // if receive and get have the same x-value, then upper will point to that value
        // otherwise, upper will point to the next highest x-value on the get array
        map<double, double>::iterator upper;
        upper = get->lower_bound(recv->first);

        // if receive and get have the same x-value, then set the coefficient as 0
        if (upper->first == recv->first) {
            coeff[recv->first] = 0;
            ++recv;
            continue;
        }
        
        // if get does not have the same x-value as receive, then calculate the interval and insert it into the map
        else if (upper != get->end()) {

            // lower points to the next x-value in get that's smaller than receive's x-value
            // upper points to the next x-value in get that's bigger than receive's x-value
            map<double, double>::iterator lower;
            lower = upper;
            --lower;

            // calculates the coefficient (with respect to the smaller x-value)
            double remainder = recv->first - lower->first;
            double interval = upper->first - lower->first;
            coeff[recv->first] = remainder / interval;
            ++recv;
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    return coeff;
}



// ONLY WORKS FOR INTERPOLATION, NOT EXTRAPOLATION

/// @brief retrieves data from another grid (MPI message passing)
/// @param receive the map that contains information on what processor has what data for the grid that wants the data
/// @param send the map that contains information on what processor has what data for the grid that needs to send the data
/// @param grid the grid that needs to send the data
/// @param coeff the coefficients 
/// @param iProc the rank of the processor that's executing the function
/// @param xPos the x-position that we are interested in getting data to and from
/// @param answer the data will be stored in this variable on the processor that wanted the data
/// @return the rank of the processor that has received the data
int getValue(map<double, double> *receive, map<double, double> *send, map<double, double> *grid, 
             map<double, double> *coefficient, int iProc, double xPos, double *answer) {

    double coeff = coefficient->at(xPos); // finds the coefficient at the requested x-value
    int rcv = findProc(receive, xPos); // rcv is the processor rank that wants to receive the data


    // if coefficient is zero, then we do not need to do any interpolating, simply just message pass
    if (coeff == 0) {

        int snd = findProc(send, xPos); // snd is the processor rank that needs to send the data

        // if the same processor contains the data for the x-value, then set answer as the data we need and return
        if (snd == rcv) {
            if (rcv == iProc) *answer = grid->at(xPos);
            return rcv;
        }

        // message passing to send the data between processors and store it in answer
        if (iProc == snd) {
            MPI_Send(&grid->at(xPos), 1, MPI_DOUBLE, rcv, 0, MPI_COMM_WORLD);
        }
        else if (iProc == rcv) {
            MPI_Recv(answer, 1, MPI_DOUBLE, snd, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        return rcv;
    }

    // this part is ran if we need to interpolate (and find weighted average)
    else {

        // upper points to x-value just above xPos
        // lower points to the x-value just below xPos
        map<double, double>::iterator upper, lower;
        upper = send->upper_bound(xPos);
        lower = upper;
        --lower;

        int sndLower = findProc(send, lower->first); // processor that contains the data at lower
        int sndUpper = findProc(send, upper->first); // processor that ocntains the data at upper
        double upperValue, lowerValue;


        // if the receive processor is the same processor that contains the data at lower, simply get that data
        if (sndLower == rcv) {
            if (rcv == iProc) lowerValue = grid->at(lower->first);
        }

        // message passing to obtain lowerValue on the rcv processor
        else {
            if (iProc == sndLower) {
                MPI_Send(&grid->at(lower->first), 1, MPI_DOUBLE, rcv, 1, MPI_COMM_WORLD);
            }
            else if (iProc == rcv) {
                MPI_Recv(&lowerValue, 1, MPI_DOUBLE, sndLower, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }

        // repeat the same process to get upperValue
        if (sndUpper == rcv) {
            if (rcv == iProc) upperValue = grid->at(upper->first);
        }
        else {
            if (iProc == sndUpper) {
                MPI_Send(&grid->at(upper->first), 1, MPI_DOUBLE, rcv, 2, MPI_COMM_WORLD);
            }
            else if (iProc == rcv) {
                MPI_Recv(&upperValue, 1, MPI_DOUBLE, sndUpper, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }

        
        // calculates the weighted average and stores it in answer
        MPI_Barrier(MPI_COMM_WORLD);
        if (iProc == rcv) {
            *answer = (coeff * lowerValue) + ((1 - coeff) * upperValue);
        }

        return rcv;
    }
}
