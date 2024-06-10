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

/*std::map<int, int> divide_list(const std::map<int, int> points, int nProcs) {
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
}*/

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

main(int argc, char **argv) {
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

    int num_points = 100;
    int seed = 42;
    std::map<int, int> points = generate_points(num_points, cummulative_start, cummulative_end, seed);

    std::map<int, int> local;
    auto iterator = points.begin();
    for (int i = 0; i < points.size(); i++) {
        if (iterator->second >= local_start && iterator->second <= local_end) {
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

    // ************************* AI CODE *************************

    int local_count = local.size();
    MPI_Bcast(&local_count, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (iProc == 0) {
        std::cout << "Processor 0 broadcasting count of local points: " << local_count << std::endl;
    }
    else {
        std::cout << "Processor " << iProc << " receiving count of local points: " << local_count << std::endl;
    }

    MPI_Finalize();
    return 0;
}