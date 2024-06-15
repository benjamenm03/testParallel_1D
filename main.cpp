// Author: Benjamen Miller, University of Michigan - Ann Arbor
// Date: 06/13/2024
// main.cpp file for testing one-dimensional examples of MPI programs for Aether

#include "main.h"
#include "mpi.h"

// COMPILE COMMAND: mpic++ -std=c++11 main.cpp functions.cpp -o main
// RUN COMMAND: mpirun -np 4 ./main

int main(int argc, char **argv) {
    // ********** Debug boolean to print testing statements **********
    bool debug_general = false;
    bool debug_print_transfer = true;
    // ***************************************************************

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

    int grid_copy_nPts = 120; // Number of points on grid_copy
    int grid_copy_start_index = 20; // Start index of grid_copy
    int grid_copy_end_index = 79; // End index of grid_copy
    std::map<double, double> grid_copy; // Initialize grid_copy map

    // Each processor generates its own chunk of grid_copy
    std::map<double, double> grid_copy_ownership = gen_grid_ref(iProc, nProcs, grid_copy, grid_copy_start_index, grid_copy_end_index, grid_copy_nPts);

    std::vector<double> packed_grid_copy = pack_map(grid_copy); // Pack grid_copy into a vector of doubles
    std::vector<double> packed_global_grid_copy(packed_grid_copy.size()); // Initialize packed_global_grid_copy vector

    // Reduce packed_grid copy from all processors into packed_global_grid_copy utilizing MPI_MAX functions
    MPI_Allreduce(&packed_grid_copy[0], &packed_global_grid_copy[0], packed_grid_copy.size(), MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    std::map<double, double> global_grid_copy = unpack_vector(packed_global_grid_copy); // Unpack packed_global_grid_copy back into a map of doubles

    std::vector<double> packed_grid_copy_ownership = pack_map(grid_copy_ownership); // Pack grid_copy_ownership into a vector of doubles
    std::vector<double> packed_global_grid_copy_ownership(packed_grid_copy_ownership.size()); // Initialize packed_global_grid_copy_ownership vector

    // Reduce packed_grid_copy_ownership from all processors into packed_global_grid_copy_ownership utilizing MPI_MAX functions
    MPI_Allreduce(&packed_grid_copy_ownership[0], &packed_global_grid_copy_ownership[0], packed_grid_copy_ownership.size(), MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    std::map<double, double> global_grid_copy_ownership = unpack_vector(packed_global_grid_copy_ownership); // Unpack packed_global_grid_copy_ownership back into a map of doubles

    if (debug_general == true) {
        // Print statements to read out the following data (print_data works with vector<double> and map<double, double>):
        print_data(iProc, global_grid_copy, "Global Grid Copy:", 0);
        print_data(iProc, global_grid_copy_ownership, "Global Grid Copy Ownership:", 0);
        print_data(iProc, temp_ref, "Local Temp Ref:", 0);
        
        // Testing get_owner() on a range of values with grid_copy
        std::map<double, double> owner_map = get_owner(iProc, grid_copy, 30, 60);
        print_data(iProc, owner_map, "grid_copy owner map for indices 30 through 60:", 0);

        // Testing get_owner() on a single index with temp_ref
        int temporary_owner = get_owner(iProc, temp_ref, 60);
        if (iProc == 0) {
            std::cout << "\nOwner of index 60 on temp_ref: " << temporary_owner << std::endl;
        }
    }

    MPI_Barrier(MPI_COMM_WORLD); // Synchronize all processors

    double index = 90; // Index to transfer
    int source_owner = get_owner(iProc, grid_copy, index);
    int dest_owner = get_owner(iProc, temp_ref, index);

    if(debug_print_transfer == true) {
        print_data(iProc, grid_copy, "Local Grid Copy:", source_owner);
        MPI_Barrier(MPI_COMM_WORLD); // Synchronize all processors
        print_data(iProc, temp_ref, "Local Temp Ref:", dest_owner);
    }

    MPI_Barrier(MPI_COMM_WORLD); // Synchronize all processors

    transfer_data(iProc, grid_copy, temp_ref, index);

    MPI_Barrier(MPI_COMM_WORLD); // Synchronize all processors

    if (debug_print_transfer == true) {
        print_data(iProc, grid_copy, "Local Grid Copy:", source_owner);
        MPI_Barrier(MPI_COMM_WORLD); // Synchronize all processors
        print_data(iProc, temp_ref, "Local Temp Ref:", dest_owner);
    }

    MPI_Finalize(); // Finalize MPI
}