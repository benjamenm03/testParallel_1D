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

std::map<int, int> gen_local_points(int iProc) {

}

int main(int argc, char **argv) {

}