// main.cpp file for testing one-dimensional examples of MPI programs for Aether

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <map>
#include <iostream>

// COMPILE COMMAND: mpic++ -std=c++11 main.cpp -o main
// RUN COMMAND: mpirun -np 4 ./main

// NEEDS ATTENTION:
// - Implement "meaningful" data exchange using message passing

bool gen_local_points(int iProc, int nProcs, std::map<int, int> &local_map, int range_end) {
    bool did_work = false; // Debugger
    int seed = 42;
    int total_length = local_map.size(); // total number of elements in the map (num_points)
    int sub_size = total_length / nProcs; // number of elements in each submap (processor num_points)
    int local_start = iProc * sub_size; // local map start index
    int local_end = local_start + sub_size - 1; // local map end index

    if (total_length <= 0) {
        std::cout << "Error: Length of local_map must be greater than 0" << std::endl;
        return did_work;
    }

    srand(seed);

    for (int i = local_start; i <= local_end; i++) {
        local_map[i] = rand() % (range_end + 1);
        did_work = true;
    }

    return did_work;
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

    int range_end = 99;
    int num_points = 100;
    std::map<int, int> local_map;
    for (int i = 0; i < num_points; i++) {
        local_map[i] = -1;
    }
    gen_local_points(iProc, nProcs, local_map, range_end);

    MPI_Barrier(MPI_COMM_WORLD);

    if (iProc == 0) {
        // Processor 0 prints first
        std::cout << "Processor " << iProc << " has the following points: " << std::endl;
        for (const auto &i : local_map) {
            std::cout << "{" << i.first << "," << i.second << "} ";
        }
        std::cout << "\nProcessor " << iProc << " has: " << local_map.size() << " points\n\n";
    }

    // Synchronize after processor 0 prints
    MPI_Barrier(MPI_COMM_WORLD);

    for (int i = 1; i < nProcs; i++) {
        if (iProc == i) {
            // Each processor waits for its turn
            std::cout << "Processor " << iProc << " has the following points: " << std::endl;
            for (const auto &j : local_map) {
                std::cout << "{" << j.first << "," << j.second << "} ";
            }
            std::cout << "\nProcessor " << iProc << " has: " << local_map.size() << " points\n\n";
        }
        // Synchronize after each processor prints
        MPI_Barrier(MPI_COMM_WORLD);
    }

    // Serialize the data from local_map.second into an array
    std::vector<int> local_array(local_map.size());
    int index = 0;
    for (std::map<int, int>::iterator i = local_map.begin(); i != local_map.end(); i++) {
        local_array[index] = i->second;
        index++;
    }

    std::vector<int> global_array(num_points);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Allreduce(&local_array[0], &global_array[0], num_points, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

    std::map<int, int> global_map;
    for (int i = 0; i < num_points; i++) {
        global_map[i] = global_array[i];
    }

    if (iProc == 0) {
        std::cout << "Global array contains the following points:" << std::endl;
        for (const auto &i : global_map) {
            std::cout << "{" << i.first << "," << i.second << "} ";
        }
    }

    MPI_Finalize();
    return 0;
}