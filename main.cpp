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
// - Need to implement "meaningful" data exchange using message passing

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

    if (iProc == 0) {
        std::cout << "\nFormat: {xLoc, val}" << std::endl;
    }

    MPI_Barrier(MPI_COMM_WORLD);

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
        if (iterator->first >= local_start && iterator->first <= local_end) {
            local[iterator->first] = iterator->second;
        }
        iterator++;
    }

    // Ensure all processors have reached this point
    MPI_Barrier(MPI_COMM_WORLD);

    if (iProc == 0) {
        // Processor 0 prints first
        std::cout << "Processor " << iProc << " has the following points: " << std::endl;
        for (const auto &i : local) {
            std::cout << "{" << i.first << "," << i.second << "} ";
        }
        std::cout << "\nProcessor " << iProc << " has: " << local.size() << " points\n\n";
    }

    // Synchronize after processor 0 prints
    MPI_Barrier(MPI_COMM_WORLD);

    for (int i = 1; i < nProcs; i++) {
        if (iProc == i) {
            // Each processor waits for its turn
            std::cout << "Processor " << iProc << " has the following points: " << std::endl;
            for (const auto &j : local) {
                std::cout << "{" << j.first << "," << j.second << "} ";
            }
            std::cout << "\nProcessor " << iProc << " has: " << local.size() << " points\n\n";
        }
        // Synchronize after each processor prints
        MPI_Barrier(MPI_COMM_WORLD);
    }

    MPI_Finalize();
    return 0;
}