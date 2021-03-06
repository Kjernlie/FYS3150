# Benchamarks

To run the benchmarks you have install MPI on your computer. How to do this is explained here: https://github.com/CompPhysics/ComputationalPhysics/tree/master/doc/Programs/ParallelizationMPI.

Once you have installed MPI, you have to be in the directory containing main.cpp, system.cpp and system.h. First, the program has to be compiled. Open the terminal, position yourself in the repository with the mentioned files, and write:

mpic++ -O3 -std=c++11 -o main.x main.cpp system.cpp -larmadillo

Now, you can run the program with 

mpirun -n 2 ./main.x benchmarks1.dat 2 250000 2.3 2.4 0.1

Here, the different inputs have the following meaning:

mpirun -n (Number of nodes to run on) ./(executable) (file, which you want to write the values to) (lattice size) (number of Monte Carlo samples to run on one node) (Initial temperature) (Final temperature) (temperature step) 


If you open benchmarks1.dat, it should contain numbers simular to 

2.4000000     -1.6457320     0.41072652     0.12418800     0.88241300    0.088973000      9.4631389     -6.5829280


Accordingly, you could write

mpirun -n 2 ./main.x benchmarks2.dat 20 250000 2.3 2.4 0.1

which should give values similar to 

2.4000000     -1.2380266      1.3995641      8.9432533     0.45400643     0.13450894      3224.5957     -495.21064

in benchmarks2.dat.





