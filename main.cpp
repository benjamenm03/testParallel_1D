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

bool gen_local_points(int iProc, int nProcs, std::map<int, int> &local_map, std::map<int, int> &ownership_map, int range_end) {
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
        ownership_map[i] = iProc;
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
        std::cout << "\nFormat: {xLoc, val}\n" << std::endl;
    }

    MPI_Barrier(MPI_COMM_WORLD);

    int range_end = 99;
    int num_points = 100;
    std::map<int, int> local_map;
    std::map<int, int> ownership_map;
    for (int i = 0; i < num_points; i++) {
        local_map[i] = -1;
        ownership_map[i] = -1;
    }

    int range_end2 = 99;
    int num_points2 = 152;
    std::map<int, int> local_map2;
    std::map<int, int> ownership_map2;
    for (int i = 0; i < num_points2; i++) {
        local_map2[i] = -1;
        ownership_map2[i] = -1;
    }

    gen_local_points(iProc, nProcs, local_map, ownership_map, range_end);
    gen_local_points(iProc, nProcs, local_map2, ownership_map2, range_end2);

    MPI_Barrier(MPI_COMM_WORLD);

    for (int i = 0; i < nProcs; i++) {
        if (iProc == i) {
            std::cout << "Processor " << iProc << " has the following points on map 1:" << std::endl;
            for (const auto &j : local_map) {
                std::cout << "{" << j.first << "," << j.second << "} ";
            }
            std::cout << "\nProcessor " << iProc << ", map 1 has: " << local_map.size() << " points\n\n";
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    for (int i = 0; i < nProcs; i++) {
        if (iProc == i) {
            std::cout << "Processor " << iProc << " has the following points on map 2:" << std::endl;
            for (const auto &j : local_map2) {
                std::cout << "{" << j.first << "," << j.second << "} ";
            }
            std::cout << "\nProcessor " << iProc << ", map 2 has: " << local_map2.size() << " points\n\n";
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    std::vector<int> local_array(local_map.size());
    int index = 0;
    for (std::map<int, int>::iterator i = local_map.begin(); i != local_map.end(); i++) {
        local_array[index] = i->second;
        index++;
    }

    std::vector<int> ownership_array(ownership_map.size());
    index = 0;
    for (std::map<int, int>::iterator i = ownership_map.begin(); i != ownership_map.end(); i++) {
        ownership_array[index] = i->second;
        index++;
    }

    std::vector<int> local_array2(local_map2.size());
    int index2 = 0;
    for (std::map<int, int>::iterator i = local_map2.begin(); i != local_map2.end(); i++) {
        local_array2[index2] = i->second;
        index2++;
    }

    std::vector<int> ownership_array2(ownership_map2.size());
    index2 = 0;
    for (std::map<int, int>::iterator i = ownership_map2.begin(); i != ownership_map2.end(); i++) {
        ownership_array2[index2] = i->second;
        index2++;
    }

    std::vector<int> global_array(num_points);
    std::vector<int> global_ownership_array(num_points);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Allreduce(&local_array[0], &global_array[0], num_points, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&ownership_array[0], &global_ownership_array[0], num_points, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

    std::vector<int> global_array2(num_points2);
    std::vector<int> global_ownership_array2(num_points2);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Allreduce(&local_array2[0], &global_array2[0], num_points2, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&ownership_array2[0], &global_ownership_array2[0], num_points2, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

    std::map<int, int> global_map;
    for (int i = 0; i < num_points; i++) {
        global_map[i] = global_array[i];
    }

    std::map<int, int> global_ownership_map;
    for (int i = 0; i < num_points; i++) {
        global_ownership_map[i] = global_ownership_array[i];
    }

    std::map<int, int> global_map2;
    for (int i = 0; i < num_points2; i++) {
        global_map2[i] = global_array2[i];
    }

    std::map<int, int> global_ownership_map2;
    for (int i = 0; i < num_points2; i++) {
        global_ownership_map2[i] = global_ownership_array2[i];
    }

    if (iProc == 0) {
        std::cout << "Global map 1 contains the following points:" << std::endl;
        for (const auto &i : global_map) {
            std::cout << "{" << i.first << "," << i.second << "} ";
        }
        std::cout << "\n\nGlobal ownership map 1:" << std::endl;
        for (const auto &i : global_ownership_map) {
            std::cout << "{" << i.first << "," << i.second << "} ";
        }
    }

    if (iProc == 0) {
        std::cout << "\n\nGlobal map 2 contains the following points:" << std::endl;
        for (const auto &i : global_map2) {
            std::cout << "{" << i.first << "," << i.second << "} ";
        }
        std::cout << "\n\nGlobal ownership map 2:" << std::endl;
        for (const auto &i : global_ownership_map2) {
            std::cout << "{" << i.first << "," << i.second << "} ";
        }
    }

    MPI_Finalize();
    return 0;
}