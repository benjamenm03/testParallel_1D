# testParallel_1D
Test repository for experimenting with one-dimensional MPI concepts before applying them to the three-dimensional Aether model.

You may need to change your include path in vscode settings, especially if you aren't working on a Mac with Homebrew.

Simple Compile Command: mpic++ -std=c++11 main.cpp functions.cpp -o main
Run Command: mpirun -np 4 ./main