// Author: Benjamen Miller, University of Michigan - Ann Arbor
// Date: 06/14/2024
// functions.cpp helper function file for testing one-dimensional examples of MPI programs for Aether

#include "main.h"

// GEN_GRID_REF:
// Modifies grid_ref (map of doubles) to be a simple single-dimensional grid of points with a logically increasing format
// within a specified range. Each processor is assigned a chunk of grid_ref to own. Returns a map of ownership, but this
// IS NOT currently being used since get_owner() is utilized instead, which is more effective.
std::map<double, double> gen_grid_ref(int iProc, int nProcs, std::map<double, double> &grid_ref, double start_index, double end_index, int nPts) {
    std::map<double, double> local_ownership_map; // Local ownership map
    int total_length = nPts; // Total number of points in grid_ref
    int sub_size = total_length / nProcs; // Number of points in each submap (processor num_points)
    int local_start = iProc * sub_size; // Local map start index
    double local_end = local_start + sub_size - 1; // Local map end index
    double step_size = (end_index - start_index + 1) / total_length; // Step size for each index

    // Fill in every index of the map with a default value of -1
    // Fill processor's chunk of grid_ref with index values to use as references for temp_ref
    // Update processor's local_ownership_map with processor ownership over its range
    for (int i = 0; i < total_length; ++i) {
        double index = start_index + i * step_size;
        grid_ref[index] = -1;
        local_ownership_map[index] = -1;
        if (i >= local_start && i <= local_end) {
            grid_ref[index] = index;
            local_ownership_map[index] = iProc;
        }
    }

    // Return local_ownership_map as function output labelling processor ownership
    return local_ownership_map;
}

// GEN_TEMP_REF:
// Modifies temp_ref (map of doubles) to be a simple single-dimensional grid of points with a sinusoidal temperature
// based on the grid_ref. Each processor is assigned a chunk of temp_ref to own. Returns a map of ownership, but this
// IS NOT currently being used since get_owner() is utilized instead, which is more effective.
std::map<double, double> gen_temp_ref(int iProc, int nProcs, std::map<double, double> &grid_ref, std::map<double, double> &temp_ref) {
    std::map<double, double> local_ownership_map; // Local ownership map
    int total_length = grid_ref.size(); // Total number of points in temp_ref
    double start_index = grid_ref.begin()->first; // Start index of grid_ref
    double end_index = grid_ref.rbegin()->first; // End index of grid_ref
    int sub_size = total_length / nProcs; // Number of points in each submap (processor num_points)
    int local_start = iProc * sub_size; // Local map start index
    double local_end = local_start + sub_size - 1; // Local map end index
    double step_size = (end_index - start_index + 1) / total_length; // Step size for each index

    // Fill in every index of the map with a default value of -1
    for (double i = start_index; i <= end_index; i += step_size) {
        temp_ref[i] = -1;
        local_ownership_map[i] = -1;
        if (i >= local_start && i <= local_end) {
            temp_ref[i] = 200 + 100 * sin(grid_ref[i] * 2 * M_PI / 100);
            local_ownership_map[i] = iProc;
        }
    }

    // Return did_work as function output boolean
    return local_ownership_map;
}

// PACK_MAP (entire map):
// Packs a map of doubles into a vector of doubles to be used in MPI handling
std::vector<double> pack_map(std::map<double, double> &map) {
    std::vector<double> packed_map;
    for (auto const &pair : map) {
        packed_map.push_back(pair.first);
        packed_map.push_back(pair.second);
    }
    return packed_map;
}

// PACK_MAP (single index):
// Packs a single index from a map of doubles into a vector of doubles to be used in MPI handling
std::vector<double> pack_map(std::map<double, double> &map, double index) {
    std::vector<double> packed_map;
    packed_map.push_back(index);
    packed_map.push_back(map[index]);
    return packed_map;
}

// PACK_MAP (range of values):
// Packs a range of values from a map of doubles into a vector of doubles to be used in MPI handling
std::vector<double> pack_map(std::map<double, double> &map, double start_index, double end_index) {
    std::vector<double> packed_map;
    for (auto const &pair : map) {
        if (pair.first >= start_index && pair.first <= end_index) {
            packed_map.push_back(pair.first);
            packed_map.push_back(pair.second);
        }
    }
    return packed_map;
}

// UNPACK_VECTOR:
// Unpacks a vector of doubles into a map of doubles to be used in MPI handling
std::map<double, double> unpack_vector(std::vector<double> &packed_map) {
    std::map<double, double> map;
    for (int i = 0; i < packed_map.size(); i += 2) {
        map[packed_map[i]] = packed_map[i + 1];
    }
    return map;
}

// GET_OWNER (single index):
// Determines the owner of a given index in a map of doubles
int get_owner(int iProc, std::map<double, double> &map, double index) {
    std::vector<int> owner_vector;
    double value = -1;
    if (map.count(index) > 0) {
        value = map[index];
    }
    if (value != -1) {
        owner_vector.push_back(iProc);
    } else {
        owner_vector.push_back(-1);
    }
    std::vector<int> global_owner_vector(owner_vector.size());
    MPI_Allreduce(&owner_vector[0], &global_owner_vector[0], owner_vector.size(), MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    int max_owner = *std::max_element(global_owner_vector.begin(), global_owner_vector.end());
    return max_owner;
}

// GET_OWNER (range of values):
// Determines the owner of a range of values in a map of doubles
std::map<double, double> get_owner(int iProc, std::map<double, double> &map, double start_index, double end_index) {
    std::map<double, double> owner_map;
    for (auto const &pair : map) {
        if (pair.first >= start_index && pair.first <= end_index) {
            if (pair.second == -1) {
                owner_map[pair.first] = -1;
            } else {
                owner_map[pair.first] = iProc;
            }
        }
    }
    std::map<double, double> global_owner_map;
    std::vector<double> local_keys;
    std::vector<double> local_values;
    for (const auto &pair : owner_map) {
        local_keys.push_back(pair.first);
        local_values.push_back(pair.second);
    }
    std::vector<double> global_keys(local_keys.size());
    std::vector<double> global_values(local_values.size());
    MPI_Allreduce(&local_keys[0], &global_keys[0], local_keys.size(), MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&local_values[0], &global_values[0], local_values.size(), MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    for (int i = 0; i < global_keys.size(); ++i) {
        global_owner_map[global_keys[i]] = global_values[i];
    }
    return global_owner_map;
}

std::pair<double, double> find_nearest_indices(std::map<double, double> &map, double index) {
    double nearest_below = -1;
    double nearest_above = -1;
    for (auto const &pair : map) {
        if (pair.first < index) {
            nearest_below = pair.first;
        } else if (pair.first > index) {
            nearest_above = pair.first;
            break;
        }
    }

    std::pair<double, double> nearest_indices = std::make_pair(nearest_below, nearest_above);
    return nearest_indices;
}

// PRINT_DATA (map<double, double>):
// Prints out a map of doubles with a header
// WARNING: This is intended to be used for "local" maps. Use this to read out data stored on a processor.
// If you want to read out data that is stored across multiple processors, use MPI_Allreduce() first.
void print_data(int iProc, std::map<double, double> &data, std::string header, int which_proc) {
    if (iProc == which_proc) {
        std::cout << "\n" << header << std::endl;
        for (const auto& pair : data) {
            std::cout << pair.first << ": " << pair.second << std::endl;
        }
    }
}

// PRINT_DATA (vector<double>):
// Prints out a vector of doubles with a header in the format of "index: value".
// WARNING: This should only be used for printing vectors that are packed and have original
// format of a map, or are intended to be seen as a map.
void print_data(int iProc, std::vector<double> &data, std::string header) {
    if (iProc == 0) {
        std::cout << "\n" << header << std::endl;
        for (int i = 0; i < data.size(); i += 2) {
            std::cout << data[i] << ": " << data[i + 1] << std::endl;
        }
    }
}

// PRINT_VECTOR:
// Prints a standard vector in csv format
void print_vector(int iProc, std::vector<double> &data) {
    if (iProc == 0) {
        std::cout << "Print Vector: " << std::endl;
        for (int i = 0; i < data.size(); i++) {
            std::cout << data[i] << ", ";
        }
    }
}

// TRANSFER_DATA (single index):
// Transfers a single index from one map to another
// Handles all cases of ownership and interpolation (that I can think of). Each if-branch describes what it's catching.
// Does not handle extrapolation for points lying outside the range of the source map.
void transfer_data(int iProc, std::map<double, double> &source_map, std::map<double, double> &dest_map, double index) {
    // Both maps contain the index
    if ((source_map.find(index) != source_map.end()) && (dest_map.find(index) != dest_map.end())) {
        int source_owner = get_owner(iProc, source_map, index);
        int dest_owner = get_owner(iProc, dest_map, index);

        // Source and destination are owned by current processor
        // Copy data from source to destination without MPI
        if (source_owner == iProc && dest_owner == iProc) {
            dest_map[index] = source_map[index];
            std::cout << "Processor " << iProc << " owns both the source and destination index. Copying data..." << std::endl;

        // Source is owned by current processor, destination is owned by another processor
        // Use MPI_Send to move data
        } else if (source_owner == iProc && dest_owner != iProc) {
            std::vector<double> packed_map(1);
            packed_map[0] = source_map[index];
            MPI_Send(&packed_map[0], 1, MPI_DOUBLE, dest_owner, 0, MPI_COMM_WORLD);
            std::cout << "Processor " << iProc << " sent data to processor " << dest_owner << std::endl;

        // Source is owned by another processor, destination is owned by current processor
        // Use MPI_Recv to receive data
        } else if (source_owner != iProc && dest_owner == iProc) {
            std::vector<double> packed_map(1);
            MPI_Recv(&packed_map[0], 1, MPI_DOUBLE, source_owner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            dest_map[index] = packed_map[0];
            std::cout << "Processor " << iProc << " received data from processor " << source_owner << std::endl;

        // Source and destination are owned by different processors
        // The current processor doesn't need to do anything!
        } else if (source_owner != iProc && dest_owner != iProc) {
            std::cout << "Processor " << iProc << " does not own the source or destination index. Skipping..." << std::endl;

        }

    // The index being transferred exists on the destination map, but not on the source map.
    // This needs interpolation!!
    } else if ((source_map.find(index) == source_map.end()) && (dest_map.find(index) != dest_map.end())) {

        // Calculates the nearest indices to the index being transferred on the source map
        std::pair<double, double> nearest_indices = find_nearest_indices(source_map, index);
        double lower_index = nearest_indices.first;
        double upper_index = nearest_indices.second;
        
        // Finds the owners of these nearest indices
        int source_owner_lower = get_owner(iProc, source_map, lower_index);
        int source_owner_upper = get_owner(iProc, source_map, upper_index);

        // Determines the owner of the destination index
        int dest_owner = get_owner(iProc, dest_map, index);

        // Calculates the interpolation coefficient
        double dist_lower = abs(lower_index - index);
        double dist_range = abs(upper_index - lower_index);
        double interpolation_coef = dist_lower / dist_range;

        // Current processor owns everything (source range and destination index)
        // Simple interpolation can be done without the use of MPI
        if (source_owner_lower == iProc && source_owner_upper == iProc && dest_owner == iProc) {
            dest_map[index] = source_map[lower_index] + interpolation_coef * (source_map[upper_index] - source_map[lower_index]);
            std::cout << "Processor " << iProc << " owns both the destination and the source range. Interpolating data..." << std::endl;

        // Current processor owns the source range but not the destination index
        // Uses MPI_Send to transfer data for destination processor to perform interpolation
        } else if (source_owner_lower == iProc && source_owner_upper == iProc && dest_owner != iProc) {
            std::vector<double> packed_data(2);
            packed_data[0] = source_map[lower_index];
            packed_data[1] = source_map[upper_index];
            MPI_Send(&packed_data[0], 2, MPI_DOUBLE, dest_owner, 0, MPI_COMM_WORLD);
            std::cout << "Processor " << iProc << " owns the source range. Sending data to processor " << dest_owner << std::endl;

        // Current processor owns the source lower index and destination index
        // Uses MPI_Recv to receive missing upper index data from source owner, and performs interpolation
        } else if (source_owner_lower == iProc && source_owner_upper != iProc && dest_owner == iProc) {
            std::vector<double> packed_data(1);
            MPI_Recv(&packed_data[0], 1, MPI_DOUBLE, source_owner_upper, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            dest_map[index] = source_map[lower_index] + interpolation_coef * (packed_data[0] - source_map[lower_index]);
            std::cout << "Processor " << iProc << " owns the destination index and lower source index. Receiving data from processor " << source_owner_upper << " and interpolating..." << std::endl;

        // Current processor owns the source upper index and destination index
        // Uses MPI_Recv to receive missing lower index data from source owner, and performs interpolation
        } else if (source_owner_lower != iProc && source_owner_upper == iProc && dest_owner == iProc) {
            std::vector<double> packed_data(1);
            MPI_Recv(&packed_data[0], 1, MPI_DOUBLE, source_owner_lower, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            dest_map[index] = packed_data[0] + interpolation_coef * (source_map[upper_index] - packed_data[0]);
            std::cout << "Processor " << iProc << " owns the destination index and upper source index. Receiving data from processor " << source_owner_lower << " and interpolating..." << std::endl;

        // Current processor owns the destination index but not the source range
        // Uses MPI_Recv to receive missing source data from source owner, and performs interpolation
        } else if (source_owner_lower != iProc && source_owner_upper != iProc && dest_owner == iProc) {
            std::vector<double> packed_data(2);
            MPI_Recv(&packed_data[0], 2, MPI_DOUBLE, source_owner_lower, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            dest_map[index] = packed_data[0] + interpolation_coef * (packed_data[1] - packed_data[0]);
            std::cout << "Processor " << iProc << " owns the destination index. Receiving data from processors " << source_owner_lower << " and " << source_owner_upper << " and interpolating..." << std::endl;

        // Current processor owns the lower source index, but not the upper source index or destination index
        // Uses MPI_Send to send lower index data to destination owner for it to perform interpolation
        } else if (source_owner_lower == iProc && source_owner_upper != iProc && dest_owner != iProc) {
            std::vector<double> packed_data(1);
            packed_data[0] = source_map[lower_index];
            MPI_Send(&packed_data[0], 1, MPI_DOUBLE, dest_owner, 0, MPI_COMM_WORLD);
            std::cout << "Processor " << iProc << " owns the lower source index. Sending data to processor " << dest_owner << std::endl;

        // Current processor owns the upper source index, but not the lower source index or destination index
        // Uses MPI_Send to send upper index data to destination owner for it to perform interpolation
        } else if (source_owner_lower != iProc && source_owner_upper == iProc && dest_owner != iProc) {
            std::vector<double> packed_data(1);
            packed_data[0] = source_map[upper_index];
            MPI_Send(&packed_data[0], 1, MPI_DOUBLE, dest_owner, 0, MPI_COMM_WORLD);
            std::cout << "Processor " << iProc << " owns the upper source index. Sending data to processor " << dest_owner << std::endl;

        // Current processor does not own the source range or destination index
        // This processor doesn't need to do anything!
        } else if (source_owner_lower != iProc && source_owner_upper != iProc && dest_owner != iProc) {
            std::cout << "Processor " << iProc << " does not own the source or destination index. Skipping..." << std::endl;

        }

        MPI_Barrier(MPI_COMM_WORLD);

    } else if ((source_map.find(index) != source_map.end()) && (dest_map.find(index) == dest_map.end())) {
        if (iProc == 0) {
            std::cout << "Error: Destination map does not contain index: " << index << std::endl;
        }
    // ***********************************************************************************************************
    } else {
        if (iProc == 0) {
            std::cout << "Error: Neither source nor destination map contain index: " << index << std::endl;
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
}

// *********** COULD IMPLEMENT A TRANSFER_DATA (RANGE OF VALUES) FUNCTION TO TRANSFER MULTIPLE INDICES ***********
// ...helper functions to do this are already implemented, just need to overload transfer_data() for ranges