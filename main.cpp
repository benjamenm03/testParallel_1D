// main.cpp file for testing one-dimensional examples of MPI programs for Aether

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <map>
#include <iostream>

// ************************* VERY MUCH INCOMPLETE CODE *************************
// COMPILE COMMAND: mpic++ -std=c++11 main.cpp -o main
// RUN COMMAND: mpirun -np 4 ./main

// NEEDS FIXING:
// - Division of points among processors. Getting weird distribution behavior... I'm tired.
// - Need to check if this is anything close to what we want or need... I'm tired.
// - Need to implement actual custom testing scripts.
// - Need to implement "meaningful" data exchange using message passing... it's 3 am.

std::map<int, int> generate_points(int num_points, int range_start, int range_end, int seed) {
    std::map<int, int> points;

    if (num_points <= 0) {
        return points;
        std::cout << "Error: Number of points must be greater than 0" << std::endl;
    }

    srand(seed);

    for (int i = 0; i < (num_points + 1); i++) {
        points[i] = rand() % (range_end - range_start + 1) + range_start;
    }

    return points;
}

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);

    int nProcs;
    int iProc;
    MPI_Comm_size(MPI_COMM_WORLD, &nProcs);

    MPI_Comm_rank(MPI_COMM_WORLD, &iProc);

    int cummulative_start = 0;
    int cummulative_end = 99;
    int sub_size = (cummulative_end - cummulative_start + 1) / nProcs;
    int local_start = cummulative_start + iProc * sub_size;
    int local_end = local_start + sub_size - 1;

    int num_points = 101;
    int seed = 42;
    std::map<int, int> points = generate_points(num_points, cummulative_start, cummulative_end, seed);

    std::map<int, int> local;
    auto iterator = points.begin();
    for (int i = 0; i < points.size(); i++) {
        if (iterator->second >= local_start && iterator->first <= local_end) {
            local[iterator->first] = iterator->second;
        }
        iterator++;
    }

    std::cout << "Processor " << iProc << " has the following points: " << std::endl;
    iterator = local.begin();
    for (int i = 0; i < local.size(); i++) {
        std::cout << iterator->second;
        if (i < local.size() - 1) {
            std::cout << ", ";
        } else {
            std::cout << std::endl;
        }
        iterator++;
    }
    
    int local_count = local.size();

    std::cout << "Processor " << iProc << " has: " << local_count << std::endl;

    MPI_Finalize();
    return 0;
}