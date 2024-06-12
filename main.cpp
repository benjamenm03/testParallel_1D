// main.cpp file for testing one-dimensional examples of MPI programs for Aether

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <map>
#include <iostream>
#include <cmath>

// COMPILE COMMAND: mpic++ -std=c++11 main.cpp -o main
// RUN COMMAND: mpirun -np 4 ./main

bool gen_grid_ref(int iProc, int nProcs, std::vector<int> &grid) {
    bool did_work = false;
    int total_length = grid.size();
    int sub_size = total_length / nProcs;
    int local_start = iProc * sub_size;
    int local_end = local_start + sub_size - 1;

    if (total_length <= 0) {
        std::cout << "Error: Size of grid must be greater than 0" << std::endl;
        return did_work;
    }

    for (int i = local_start; i <= local_end; i++) {
        grid[i] = i;
        did_work = true;
    }

    return did_work;
}

bool gen_temp_ref(int iProc, int nProcs, std::vector<int> &grid_ref, std::vector<int> &grid) {
    bool did_work = false;
    int total_length = grid.size();
    int sub_size = total_length / nProcs;
    int local_start = iProc * sub_size;
    int local_end = local_start + sub_size - 1;

    if (total_length <= 0) {
        std::cout << "Error: Size of grid must be greater than 0" << std::endl;
        return did_work;
    }

    for (int i = local_start; i <= local_end; i++) {
        grid[i] = 200 + 100 * sin(i * 2 * M_PI / 100);
    }
}

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);

    int nProcs, iProc;
    MPI_Comm_size(MPI_COMM_WORLD, &nProcs);
    MPI_Comm_rank(MPI_COMM_WORLD, &iProc);

    MPI_Barrier(MPI_COMM_WORLD);

    int nPts = 100;
    int start_index = 0;
    int end_index = 99;
    std::vector<int> grid_ref(nPts);


}