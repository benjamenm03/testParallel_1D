#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

bool generate_points(int *array, int num_points, int range_start, int range_end, int seed) {
    bool did_work = false;

    if (num_points <= 0) {
        return did_work;
    }

    srand(seed);
    array[0] = num_points;
    for (int i = 1; i < num_points + 1; i++) {
        array[i] = rand() % (range_end - range_start + 1) + range_start;
        did_work = true;
    }
    return did_work;
}

bool divide_list(int *array, int nProcs, int *owner_array) {
    bool did_work = false;

    int array_size = array[0];
    int sub_size = array_size / nProcs;
    owner_array[0] = array_size;

    for (int i = 0; i < nProcs; i++) {
        int j = 0;
        while (j < sub_size) {
            owner_array[i * sub_size + j + 1] = i;
            j++;
        }
        if (i == (nProcs - 1)) {
            did_work = true;
        }
    }
    return did_work;
}

main() {
    
}