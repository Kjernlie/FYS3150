# Project 4


## Structure

In the Project4 directory you will find the source code of the project. How to run the project is described in the section below. 

In the benchmarks direcory you will find a couple of benchmarks of the project, together with an explanation on how to run them.

The plotting directory contains all the pyhon-files used for plotting in this project.

Lastly, you will find the project report in the report repository.

## How to run the code

To run the code you have to install MPI on your computer. How to do this is explained here: https://github.com/CompPhysics/ComputationalPhysics/tree/master/doc/Programs/ParallelizationMPI.

Once you have installed MPI, you have to be in the directory containing main.cpp, system.cpp and system.h. First, the program has to be compiled. Open the terminal, position yourself in the repository with the mentioned files, and write:

mpic++ -O3 -std=c++11 -o main.x main.cpp system.cpp -larmadillo

Now, you can run the program with 

mpirun -n 2 ./main.x test.dat 2 250000 2.3 2.4 0.1

Here, the different inputs have the following meaning:

mpirun -n (Number of nodes to run on) ./(executable) (file, which you want to write the values to) (lattice size) (number of Monte Carlo samples to run on one node) (Initial temperature) (Final temperature) (temperature step) 






