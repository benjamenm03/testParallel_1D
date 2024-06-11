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
    bool did_work = false;
    int total_length = local_map.size();
    int sub_size = total_length / nProcs;
    int local_start = iProc * sub_size;
    int local_end = local_start + sub_size - 1;

    for (int i = local_start; i <= local_end; i++) {
        local_map[i] = rand() % (range_end + 1);
    }

    return did_work;
}

int main(int argc, char **argv) {

}