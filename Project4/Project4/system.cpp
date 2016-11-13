#include "mpi.h"
#include "system.h"
#include <cmath>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <armadillo>
#include <random>
#include <string>
#include <sstream>

using namespace std;
using namespace arma;

System::System(int N_spins):
    m_N_spins(N_spins),
    m_spinMatrix(zeros<mat>(N_spins,N_spins)),
    m_accepted_states(0),
    m_tol(1)

{

}


// --------------------------------------------------------------------------------------------------



void System::RunSystem(string filename, int MC_cycles, double initial_temp, double final_temp, double temp_step, int rankProcess, int NProcesses)
{
    // Start new MC sampling by looping over T first
    for (double temp = initial_temp; temp <= final_temp; temp+=temp_step)
    {

        //m_accepted_states = 0;

        vec local_expectation_vals = zeros<mat>(5);

        // start Monte Carlo computation and get local expectation values
        MetropolisSampling(MC_cycles, temp, local_expectation_vals, rankProcess, NProcesses);


        // Find the total average
        vec total_expectation_vals = zeros<mat>(5);
        for (int i=0; i<5; i++){
            MPI_Reduce(&local_expectation_vals[i], &total_expectation_vals[i], 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        }

        if (rankProcess == 0) output(filename, m_file1, MC_cycles*NProcesses, temp, total_expectation_vals);


    }
}


// --------------------------------------------------------------------------------------------------



void System::InitiliazeLattice(double& energy, double& magneticMoment)
{
    // Initialize the seed and calle the Mersenne algorithm
    random_device rd;
    mt19937_64 gen(rd());

    // Setup the uniform distribution for x in [0,1]
    uniform_real_distribution<double> distribution(0.0,1.0);

    // setup spin matrix and initial magnetization
    for(int x = 0; x < m_N_spins; x++){
        for(int y = 0; y < m_N_spins; y++){
            m_spinMatrix(x,y) = 1.0; // Spin orientation for the ground state
            magneticMoment += (double) m_spinMatrix(x,y);

        }
    }

// Random init

//    for (int x = 0; x < m_N_spins; x++){
//        for (int y = 0; y < m_N_spins; y++){
//            m_spinMatrix(x,y) = distribution(gen) < 0.5 ? 1 : -1;
//            magneticMoment += (double) m_spinMatrix(x,y);
//        }
//    }

    // setup initial energy
    for(int x = 0; x < m_N_spins; x++){
        for(int y = 0; y < m_N_spins; y++){
            energy -= (double) m_spinMatrix(x,y)*
                    (m_spinMatrix(periodic(x,m_N_spins,-1),y) +
                     m_spinMatrix(x,periodic(y,m_N_spins,-1)));
        }

    }

//    expectationValues(0) += energy;
//    expectationValues(1) += energy*energy;
//    expectationValues(2) += magneticMoment;
//    expectationValues(3) += magneticMoment*magneticMoment;
//    expectationValues(4) += fabs(magneticMoment);

}



// --------------------------------------------------------------------------------------------------



void System::MetropolisSampling(int MC_cycles, double temp, vec &expectationValues, int rankProcess, int NProcesses)
{

    // Initialize the seed and calle the Mersenne algorithm
    random_device rd;
    mt19937_64 gen(rd());

    // Setup the uniform distribution for x in [0,1]
    uniform_real_distribution<double> distribution(0.0,1.0);


    // Initialize energy and magnetization
    double energy = 0.;
    double magneticMoment = 0.;
    m_accepted_states = 0;
//    bool burning = 0;
//    double energyOld = 1000;

    // Initialize array for expectation values
    InitiliazeLattice(energy, magneticMoment);

    // setup array for possible energy chances

    vec energyDifference = zeros<mat>(17);
    //double energyOld = calculateenergy(m_spinMatrix, m_N_spins);
    for (int de=-8; de <= 8; de+=4) energyDifference(de+8) = exp(-de/temp);

    for (int cycles = 1; cycles <= MC_cycles; cycles++){

        // The sweep over the lattice, looping over all spin sites
        for (int x=0; x < m_N_spins; x++){
            for (int y=0; y < m_N_spins; y++){

                int ix = (int) (distribution(gen)*(double)m_N_spins);
                int iy = (int) (distribution(gen)*(double)m_N_spins);



                int delta_E = 2*m_spinMatrix(ix,iy)*
                        (m_spinMatrix(ix,periodic(iy,m_N_spins,-1)) +
                         m_spinMatrix(periodic(ix,m_N_spins,-1),iy) +
                         m_spinMatrix(ix,periodic(iy,m_N_spins,1)) +
                         m_spinMatrix(periodic(ix,m_N_spins,1),iy));



//                int delta_E = calculate_deltaE(m_spinMatrix, ix, iy, m_N_spins);

//                double energyNew = calculateenergy(m_spinMatrix, m_N_spins);
//                double delta_E = energyNew-energyOld;


                if (distribution(gen) <= energyDifference(delta_E+8)) {//exp(-delta_E*1/temp)) {
                    m_spinMatrix(ix,iy) *= -1.0; // flip one spin and accept new spin config
                    magneticMoment += (double) 2*m_spinMatrix(ix,iy);
                    energy += (double) delta_E;
                    m_accepted_states++;
                    //energyOld = energyNew;
                }
            }
        }
        // update expectation value for local node
        expectationValues(0) += energy;
        expectationValues(1) += energy*energy;
        expectationValues(2) += magneticMoment;
        expectationValues(3) += magneticMoment*magneticMoment;
        expectationValues(4) += fabs(magneticMoment);


        if (cycles >= 10000 && cycles % 100 == 0){
            intermediate_output(cycles*NProcesses, temp, expectationValues, MC_cycles, rankProcess);
        }

    }
}

// --------------------------------------------------------------------------------------------------



inline int System::periodic(int i, int limit, int add){
    return (i+limit+add) % limit;
}


// --------------------------------------------------------------------------------------------------


void System::output(string filename, ofstream &file, int MC_cycles, double temp, vec expectationValues)
{
    if(!file.good()) {
        file.open(filename.c_str(), ofstream::out);
        if(!file.good()) {
            cout << "Error opening file " << filename << ". Aborting!" << endl;
            terminate();
        }
    }

    double norm = 1.0/((double) (MC_cycles));  // divided by  number of cycles
    double E_expectationValues = expectationValues(0)*norm; // energy
    double E2_expectationValues = expectationValues(1)*norm;  // energy squared
    double M_expectationValues = expectationValues(2)*norm;  // Magnetic moments
    double M2_expectationValues = expectationValues(3)*norm;   //  Magnetic moments squared
    double Mabs_expectationValues = expectationValues(4)*norm;  // Absolute value magnetic moment

    // all expectation values are per spin, divide by 1/NSpins/NSpins
    double Evariance = (E2_expectationValues- E_expectationValues*E_expectationValues)/m_N_spins/m_N_spins;
    double Mvariance = (M2_expectationValues - M_expectationValues*M_expectationValues)/m_N_spins/m_N_spins;

    double energyPrSpin = E_expectationValues/m_N_spins/m_N_spins;
    double specHeatPrSpin = Evariance/temp/temp;
    double magMomPrSpin = M_expectationValues/m_N_spins/m_N_spins;
    double susceptPrSpin = Mvariance/temp;
    double absMagMomPrSpin = Mabs_expectationValues/m_N_spins/m_N_spins;

    file << setiosflags(ios::showpoint | ios::uppercase);

    file << setw(15) << setprecision(8) << temp;
    file << setw(15) << setprecision(8) << energyPrSpin;
    file << setw(15) << setprecision(8) << specHeatPrSpin;
    //m_file1 << setw(15) << setprecision(8) << magMomPrSpin;
    file << setw(15) << setprecision(8) << susceptPrSpin;
    file << setw(15) << setprecision(8) << absMagMomPrSpin;
    file << setw(15) << setprecision(8) << m_accepted_states/((double) MC_cycles*m_N_spins*m_N_spins);
    file << "\n";

}


// --------------------------------------------------------------------------------------------------




void System::intermediate_output(int cycles, double temp, vec local_expectation_vals, int MC_cycles, int rankProcess)
{

    vec total_expectation_vals = zeros<mat>(5);
    for (int i=0; i<5; i++){
        MPI_Reduce(&local_expectation_vals[i], &total_expectation_vals[i], 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    }

    ostringstream oss;
    oss << "intermediate_" << m_N_spins << "L_" << MC_cycles << "cycles.dat";
    string filename = oss.str();

    if (rankProcess == 0) output(filename, m_file2, cycles, temp, total_expectation_vals);
}




















//int System::calculate_deltaE(int ix, int iy, int N_spins){

//    int delta_E = 2*m_spinMatrix(ix,iy)*
//            (m_spinMatrix(ix,periodic(iy,N_spins,-1)) +
//             m_spinMatrix(periodic(ix,N_spins,-1),iy) +
//             m_spinMatrix(ix,periodic(iy,N_spins,1)) +
//             m_spinMatrix(periodic(ix,N_spins,1),iy));

//    return delta_E;
//}



//int System::calculateEnergy() {
//    double energy = 0;

//    for (int i = 0; i<N_spins; i++){
//        for (int j = 0; j<N_spins; j++){
//            int iNext = i+1;
//            int iLast = i-1;
//            int jNext = j+1;
//            int jLast = j-1;
//            if (i == N_spins-1) iNext = 0;
//            if (i == 0) iLast = N_spins-1;
//            if (j == N_spins-1) jNext = 0;
//            if (j == 0) jLast = N_spins-1;

//            energy += m_spinMatrix(i,j)*(m_spinMatrix(iNext,j) +
//                                       m_spinMatrix(iLast,j) +
//                                       m_spinMatrix(i,jLast) +
//                                       m_spinMatrix(i,jNext));

//        }
//    }
//    //cout << N_spins << endl;
//    return -energy*0.5;

//}







