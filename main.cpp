// main.cpp file for testing one-dimensional examples of MPI programs for Aether

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <map>
#include <iostream>
#include <cmath>
#include <vector>

// COMPILE COMMAND: mpic++ -std=c++11 main.cpp -o main
// RUN COMMAND: mpirun -np 4 ./main

std::map<double, double> gen_grid_ref(int iProc, int nProcs, std::map<double, double> &grid_ref, double start_index, double end_index, int nPts) {
    std::map<double, double> local_ownership_map;
    int total_length = nPts; // Total number of points in grid_ref
    int sub_size = total_length / nProcs; // Number of points in each submap (processor num_points)
    int local_start = iProc * sub_size; // Local map start index
    double local_end = local_start + sub_size - 1; // Local map end index
    double step_size = (end_index - start_index + 1) / total_length; // Step size for each index

    // Fill in every index of the map with a default value of -1
    for (int i = 0; i < total_length; ++i) {
        double index = start_index + i * step_size;
        grid_ref[index] = -1;
        local_ownership_map[index] = -1;
        if (i >= local_start && i <= local_end) {
            grid_ref[index] = index;
            local_ownership_map[index] = iProc;
        }
    }

    // Return did_work as function output boolean
    return local_ownership_map;
}

std::map<double, double> gen_temp_ref(int iProc, int nProcs, std::map<double, double> &grid_ref, std::map<double, double> &temp_ref) {
    std::map<double, double> local_ownership_map;
    int total_length = grid_ref.size(); // Total number of points in temp_ref
    double start_index = grid_ref.begin()->first; // Start index of grid_ref
    double end_index = grid_ref.rbegin()->first; // End index of grid_ref
    int sub_size = total_length; // Number of points in each submap (processor num_points)
    int local_start = iProc * sub_size; // Local map start index
    double local_end = local_start + sub_size - 1; // Local map end index
    double step_size = (end_index - start_index + 1) / total_length; // Step size for each index

    for (double i = start_index; i <= end_index; i += step_size) {
        temp_ref[i] = -1;
        local_ownership_map[i] = -1;
    }

    // Generate processor's chunk of temp_ref
    for (double i = local_start; i <= local_end; i += step_size) {
        temp_ref[i] = 200 + 100 * sin(grid_ref[i] * 2 * M_PI / 100);
        local_ownership_map[i] = iProc;
    }

    // Return did_work as function output boolean
    return local_ownership_map;
}

std::vector<double> pack_map(std::map<double, double> &map) {
    std::vector<double> packed_map(map.size() * 2);
    for (auto const &pair : map) {
        packed_map.push_back(pair.first);
        packed_map.push_back(pair.second);
    }
    return packed_map;
}

std::map<double, double> unpack_vector(std::vector<double> &packed_map) {
    std::map<double, double> map;
    for (int i = 0; i < packed_map.size(); i += 2) {
        map[packed_map[i]] = packed_map[i + 1];
    }
    return map;
}

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv); // Initialize MPI

    int nProcs, iProc; // Number of processors, processor rank
    MPI_Comm_size(MPI_COMM_WORLD, &nProcs); // Get number of processors
    MPI_Comm_rank(MPI_COMM_WORLD, &iProc); // Get processor rank

    MPI_Barrier(MPI_COMM_WORLD); // Synchronize all processors

    int grid_ref_nPts = 100; // Number of points on grid_ref
    int grid_ref_start_index = 0; // Start index of grid_ref
    int grid_ref_end_index = 99; // End index of grid_ref
    std::map<double, double> grid_ref; // Initialize grid_ref map
    std::map<double, double> temp_ref; // Initialize temp_ref map

    // Each processor generates its own chunk of grid_ref and temp_ref
    std::map<double, double> grid_ref_ownership = gen_grid_ref(iProc, nProcs, grid_ref, grid_ref_start_index, grid_ref_end_index, grid_ref_nPts);
    std::map<double, double> temp_ref_ownership = gen_temp_ref(iProc, nProcs, grid_ref, temp_ref);

    int grid_copy_nPts = 120;
    int grid_copy_start_index = 20;
    int grid_copy_end_index = 79;
    std::map<double, double> grid_copy;

    std::map<double, double> grid_copy_ownership = gen_grid_ref(iProc, nProcs, grid_copy, grid_copy_start_index, grid_copy_end_index, grid_copy_nPts);

    std::vector<double> packed_grid_copy = pack_map(grid_copy);
    std::vector<double> packed_global_grid_copy(packed_grid_copy.size());
    MPI_Allreduce(&packed_grid_copy[0], &packed_global_grid_copy[0], packed_grid_copy.size(), MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    std::map<double, double> global_grid_copy = unpack_vector(packed_global_grid_copy);

    std::vector<double> packed_grid_copy_ownership = pack_map(grid_copy_ownership);
    std::vector<double> packed_global_grid_copy_ownership(packed_grid_copy_ownership.size());
    MPI_Allreduce(&packed_grid_copy_ownership[0], &packed_global_grid_copy_ownership[0], packed_grid_copy_ownership.size(), MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    std::map<double, double> global_grid_copy_ownership = unpack_vector(packed_global_grid_copy_ownership);

    if (iProc == 0) {
        std::cout << "Global Grid Copy:" << std::endl;
        for (const auto& pair : global_grid_copy) {
            std::cout << pair.first << ": " << pair.second << std::endl;
        }

        std::cout << "Global Grid Copy Ownership:" << std::endl;
        for (const auto& pair : global_grid_copy_ownership) {
            std::cout << pair.first << ": " << pair.second << std::endl;
        }
    }

    MPI_Finalize(); // Finalize MPI
}