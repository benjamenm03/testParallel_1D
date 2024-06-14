#include "main.h"

// GEN_GRID_REF:
// Modifies grid_ref (map of doubles) to be a simple single-dimensional grid of points with a logically increasing format
// within a specified range. Each processor is assigned a chunk of grid_ref to own.
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
// based on the grid_ref. Each processor is assigned a chunk of temp_ref to own.
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

// DET_OWNER (single index):
// Determines the owner of a given index in a map of doubles
int det_owner(int iProc, std::map<double, double> &map, double index) {
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

// DET_OWNER (range of values):
// Determines the owner of a range of values in a map of doubles
std::map<double, double> det_owner(int iProc, std::map<double, double> &map, double start_index, double end_index) {
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

// PRINT_DATA (map<double, double>):
// Prints out a map of doubles with a header
void print_data(int iProc, std::map<double, double> &data, std::string header) {
    if (iProc == 0) {
        std::cout << "\n" << header << std::endl;
        for (const auto& pair : data) {
            std::cout << pair.first << ": " << pair.second << std::endl;
        }
    }
}

// PRINT_DATA (vector<double>):
// Prints out a vector of doubles with a header in the format of "index: value".
// This should only be used for printing vectors that are packed and have original
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