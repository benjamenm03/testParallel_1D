// main.cpp file for testing one-dimensional examples of MPI programs for Aether

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <map>
#include <iostream>

// ************************* VERY MUCH INCOMPLETE CODE *************************

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

std::map<int, int> divide_list(const std::map<int, int> points, int nProcs) {
    std::map<int, int> owners;

    int map_size = points.size();
    int sub_size = map_size / nProcs;

    for (int i = 0; i < nProcs; i++) {
        int j = 0;
        while (j < sub_size) {
            j++;
        }
    }

    return owners;
}

main() {
    
}