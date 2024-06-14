#ifndef MAIN_H
#define MAIN_H

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <map>
#include <iostream>
#include <cmath>
#include <vector>
#include <string>

// GEN_GRID_REF:
// Modifies grid_ref (map of doubles) to be a simple single-dimensional grid of points with a logically increasing format
// within a specified range. Each processor is assigned a chunk of grid_ref to own.
std::map<double, double> gen_grid_ref(int iProc, int nProcs, std::map<double, double> &grid_ref, double start_index, double end_index, int nPts);

// GEN_TEMP_REF:
// Modifies temp_ref (map of doubles) to be a simple single-dimensional grid of points with a sinusoidal temperature
// based on the grid_ref. Each processor is assigned a chunk of temp_ref to own.
std::map<double, double> gen_temp_ref(int iProc, int nProcs, std::map<double, double> &grid_ref, std::map<double, double> &temp_ref);

// PACK_MAP (entire map):
// Packs a map of doubles into a vector of doubles to be used in MPI handling
std::vector<double> pack_map(std::map<double, double> &map);

// PACK_MAP (single index):
// Packs a single index from a map of doubles into a vector of doubles to be used in MPI handling
std::vector<double> pack_map(std::map<double, double> &map, double index);

// PACK_MAP (range of values):
// Packs a range of values from a map of doubles into a vector of doubles to be used in MPI handling
std::vector<double> pack_map(std::map<double, double> &map, double start_index, double end_index);

// UNPACK_VECTOR:
// Unpacks a vector of doubles into a map of doubles to be used in MPI handling
std::map<double, double> unpack_vector(std::vector<double> &packed_map);

// DET_OWNER (single index):
// Determines the owner of a given index in a map of doubles
int det_owner(int iProc, std::map<double, double> &map, double index);

// DET_OWNER (range of values):
// Determines the owner of a range of values in a map of doubles
std::map<double, double> det_owner(int iProc, std::map<double, double> &map, double start_index, double end_index);

// PRINT_DATA (map<double, double>):
// Prints out a map of doubles with a header
void print_data(int iProc, std::map<double, double> &data, std::string header);

// PRINT_DATA (vector<double>):
// Prints out a vector of doubles with a header in the format of "index: value".
// This should only be used for printing vectors that are packed and have original
// format of a map, or are intended to be seen as a map.
void print_data(int iProc, std::vector<double> &data, std::string header);

// PRINT_VECTOR:
// Prints a standard vector in csv format
void print_vector(int iProc, std::vector<double> &data);

#endif