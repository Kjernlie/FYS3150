#include "system.h"
#include <cmath>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <armadillo>
#include <random>
#include <string>

using namespace std;
using namespace arma;

System::System()
{

}


void System::RunSystem(int MC_cycles, int N_spins, double initial_temp, double final_temp, double temp_step)
{
    // Start new MC sampling by looping over T first
    for (double Temperature = initial_temp; Temperature <= final_temp; Temperature+=temp_step)
    {
        vec ExpectationValues = zeros<mat>(5);

        // Start MC computation
        MetrepolisSampling(N_spins,MC_cycles, Temperature, ExpectationValues);
        output("test.dat", N_spins, MC_cycles, Temperature, ExpectationValues);
    }
}

void System::InitiliazeLattice(int N_spins, mat &SpinMatrix, double& Energy, double& MagneticMoment)
{
    // setup spin matrix and initial magnetization
    for(int x = 0; x < N_spins; x++){
        for(int y = 0; y < N_spins; y++){
            SpinMatrix(x,y) = 1.0; // Spin orientation for the ground state
            MagneticMoment += (double) SpinMatrix(x,y);

        }
    }

    // setup initial energy
    for(int x = 0; x < N_spins; x++){
        for(int y = 0; y < N_spins; y++){
            Energy -= (double) SpinMatrix(x,y)*
                    (SpinMatrix(periodic(x,N_spins,-1),y) +
                     SpinMatrix(x,periodic(y,N_spins,-1)));
        }
    }
}



void System::MetrepolisSampling(int N_spins, int MC_cycles, double Temperature, vec &ExpectationValues)
{

    // Initialize the seed and calle the Mersenne algorithm
    random_device rd;
    mt19937_64 gen(rd());

    // Setup the uniform distribution for x in [0,1]
    uniform_real_distribution<double> distribution(0.0,1.0);

    // Allocate memory for spin matrix
    mat SpinMatrix = zeros<mat>(N_spins,N_spins);

    // Initialize energy and magnetization
    double Energy = 0.;
    double MagneticMoment = 0.;

    // Initialize array for expectation values
    InitiliazeLattice(N_spins, SpinMatrix, Energy, MagneticMoment);

    // setup array for possible energy chances
    vec EnergyDifference = zeros<mat>(17);
    double EnergyOld = calculateEnergy(SpinMatrix, N_spins);
    //for (int de=-8; de <= 8; de+=4) EnergyDifference(de+8) = exp(-de/Temperature);
    for (int cycles = 1; cycles <= MC_cycles; cycles++){

        // The sweep over the lattice, looping over all spin sites
        for (int x=0; x < N_spins; x++){
            for (int y=0; y < N_spins; y++){

                int ix = (int) (distribution(gen)*(double)N_spins);
                int iy = (int) (distribution(gen)*(double)N_spins);

//                int delta_E = 2*SpinMatrix(ix,iy)*
//                        (SpinMatrix(ix,periodic(iy,N_spins,-1)) +
//                         SpinMatrix(periodic(ix,N_spins,-1),iy) +
//                         SpinMatrix(ix,periodic(iy,N_spins,1)) +
//                         SpinMatrix(periodic(ix,N_spins,1),iy));

                double EnergyNew = calculateEnergy(SpinMatrix, N_spins);
                double delta_E = EnergyNew-EnergyOld;

                if (distribution(gen) <= exp(-delta_E*1/Temperature)) {//EnergyDifference(delta_E+8)) {
                    SpinMatrix(ix,iy) *= -1.0; // flip one spin and accept new spin config
                    MagneticMoment += (double) 2*SpinMatrix(ix,iy);
                    Energy += (double) delta_E;
                }
            }
        }
        // update expectation value for local node
        ExpectationValues(0) += Energy;
        ExpectationValues(1) += Energy*Energy;
        ExpectationValues(2) += MagneticMoment;
        ExpectationValues(3) += MagneticMoment*MagneticMoment;
        ExpectationValues(4) += fabs(MagneticMoment);
    }
}

inline int System::periodic(int i, int limit, int add){
    return (i+limit+add) % limit;
}

void System::output(string filename, int N_spins, int MC_cycles, double temperature, vec ExpectationValues)
{
    if(!m_file.good()) {
        m_file.open(filename.c_str(), ofstream::out);
        if(!m_file.good()) {
            cout << "Error opening file " << filename << ". Aborting!" << endl;
            terminate();
        }
    }

    double norm = 1.0/((double) (MC_cycles));  // divided by  number of cycles
    double E_ExpectationValues = ExpectationValues(0)*norm;
    double E2_ExpectationValues = ExpectationValues(1)*norm;
    double M_ExpectationValues = ExpectationValues(2)*norm;
    double M2_ExpectationValues = ExpectationValues(3)*norm;
    double Mabs_ExpectationValues = ExpectationValues(4)*norm;
    // all expectation values are per spin, divide by 1/NSpins/NSpins
    double Evariance = (E2_ExpectationValues- E_ExpectationValues*E_ExpectationValues)/N_spins/N_spins;
    double Mvariance = (M2_ExpectationValues - Mabs_ExpectationValues*Mabs_ExpectationValues)/N_spins/N_spins;
    m_file << setiosflags(ios::showpoint | ios::uppercase);
    m_file << setw(15) << setprecision(8) << temperature;
    m_file << setw(15) << setprecision(8) << E_ExpectationValues/N_spins/N_spins;
    m_file << setw(15) << setprecision(8) << Evariance/temperature/temperature;
    m_file << setw(15) << setprecision(8) << M_ExpectationValues/N_spins/N_spins;
    m_file << setw(15) << setprecision(8) << Mvariance/temperature;
    m_file << setw(15) << setprecision(8) << Mabs_ExpectationValues/N_spins/N_spins;

}

int System::calculateEnergy(mat spinMatrix, int N_spins) {
    double Energy = 0;

    for (int i = 0; i<N_spins; i++){
        for (int j = 0; j<N_spins; j++){
            int iNext = i+1;
            int iLast = i-1;
            int jNext = j+1;
            int jLast = j-1;
            if (i == N_spins-1) iNext = 0;
            if (i == 0) iLast = N_spins-1;
            if (j == N_spins-1) jNext = 0;
            if (j == 0) jLast = N_spins-1;

            Energy += spinMatrix(i,j)*(spinMatrix(iNext,j) +
                                       spinMatrix(iLast,j) +
                                       spinMatrix(i,jLast) +
                                       spinMatrix(i,jNext));

        }
    }
    //cout << N_spins << endl;
    return -Energy*0.5;

}







