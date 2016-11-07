#include "mpi.h"
#include <iostream>
#include "system.h"
#include <string>

using namespace std;

int main(int argc, char* argv[])
{
//    int N_spins = 2;
//    int MC_cycles = 1e5;
//    double initial_temp = 0.0;
//    double final_temp = 1.0;
//    double temp_step = 0.1;

    // mpi stufff -------------------------
    int NProcesses;
    int rankProcess;
    string filename;
    int N_spins;
    int MC_cycles;
    double initial_temp;
    double final_temp;
    double temp_step;


    //  MPI initializations
    MPI_Init (&argc, &argv);
    MPI_Comm_size (MPI_COMM_WORLD, &NProcesses);
    MPI_Comm_rank (MPI_COMM_WORLD, &rankProcess);
    if (rankProcess == 0 && argc <= 5) {
        cout << "Bad Usage: " << argv[0] <<
            " read output file, Number of spins, MC cycles, initial and final temperature and tempurate step" << endl;
        exit(1);
    }
    if ((rankProcess == 0) && (argc > 1)) {
        filename=argv[1];
        N_spins = atoi(argv[2]);
        MC_cycles = atoi(argv[3]);
        initial_temp = atof(argv[4]);
        final_temp = atof(argv[5]);
        temp_step = atof(argv[6]);
    }

    // broadcast to all nodes common variables since only master node reads from command line
    MPI_Bcast (&MC_cycles, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast (&N_spins, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast (&initial_temp, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast (&final_temp, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast (&temp_step, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

// ----------------------------------------

    System test;

    double  time_start, time_end, total_time;
    time_start = MPI_Wtime();


    test.RunSystem(filename, MC_cycles, N_spins, initial_temp, final_temp, temp_step, rankProcess, NProcesses);


    time_end = MPI_Wtime();
    total_time = time_end - time_start;

    if (rankProcess == 0){
        cout << "Time = " <<  total_time  << " on number of processors: "  << NProcesses  << endl;
    }

    MPI_Finalize ();
    //test.RunSystem(filename, MC_cycles, N_spins, initial_temp, final_temp, temp_step);



    return 0;
}
