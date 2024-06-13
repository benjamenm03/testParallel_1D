// main.cpp file for testing one-dimensional examples of MPI programs for Aether

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <map>
#include <iostream>
#include <cmath>

// COMPILE COMMAND: mpic++ -std=c++11 main.cpp -o main
// RUN COMMAND: mpirun -np 4 ./main

bool gen_grid_ref(int iProc, int nProcs, std::map<double, double> &grid_ref, double start_index, double end_index) {
    bool did_work = false; // Debugger
    int total_length = grid_ref.size(); // Total number of points in grid_ref
    int sub_size = total_length / nProcs; // Number of points in each submap (processor num_points)
    int local_start = iProc * sub_size; // Local map start index
    double local_end = local_start + sub_size - 1; // Local map end index
    double step_size = (end_index - start_index + 1) / total_length; // Step size for each index
    std::cout << "Step size: " << step_size << "\n"; // Print step size checker

    // Check if total_length is less than or equal to 0
    if (total_length <= 0) {
        std::cout << "Error: Size of grid_ref must be greater than 0" << std::endl;
        return did_work;
    }

    // Generate processor's chunk of grid_ref
    for (double i = local_start; i <= local_end; i += step_size) {
        grid_ref[i] = i;
        did_work = true;
    }

    // Return did_work as function output boolean
    return did_work;
}

bool gen_temp_ref(int iProc, int nProcs, std::map<double, double> &grid_ref, std::map<double, double> &temp_ref) {
    bool did_work = false; // Debugger
    int total_length = temp_ref.size(); // Total number of points in temp_ref
    double start_index = grid_ref.begin()->first; // Start index of grid_ref
    double end_index = grid_ref.rbegin()->first; // End index of grid_ref
    int sub_size = total_length; // Number of points in each submap (processor num_points)
    int local_start = iProc * sub_size; // Local map start index
    double local_end = local_start + sub_size - 1; // Local map end index
    double step_size = (end_index - start_index + 1) / total_length; // Step size for each index
    std::cout << "Step size: " << step_size << "\n"; // Print step size checker

    // Check if total_length is less than or equal to 0
    if (total_length <= 0) {
        std::cout << "Error: Size of temp_ref must be greater than 0" << std::endl;
        return did_work;
    }

    // Generate processor's chunk of temp_ref
    for (double i = local_start; i <= local_end; i += step_size) {
        temp_ref[i] = 200 + 100 * sin(grid_ref[i] * 2 * M_PI / 100);
        did_work = true;
    }

    // Return did_work as function output boolean
    return did_work;
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
    gen_grid_ref(iProc, nProcs, grid_ref, grid_ref_start_index, grid_ref_end_index);
    gen_temp_ref(iProc, nProcs, grid_ref, temp_ref);

    MPI_Finalize(); // Finalize MPI
}