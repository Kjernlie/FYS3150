#include "vmc.h"
#include <iostream>
#include <armadillo>
#include <random>

using namespace std;



VMC::VMC() :
    nDimensions(3),
    omega(1),
    stepLength(0.1),
    nParticles(2),
    alpha(1),
    nCycles(100000)
{
}


void VMC::runMCIntegration()
{
    rOld = zeros<mat>(nParticles, nDimensions);
    rNew = zeros<mat>(nParticles, nDimensions);
    double waveFuncOld = 0;
    double waveFuncNew = 0;
    double energySum = 0;
    double deltaE;

    // Initialize the seed and calle the Mersenne algorithm
    random_device rd;
    mt19937_64 gen(rd());

    // Setup the uniform distribution for x in [0,1]
    uniform_real_distribution<double> distribution(0.0,1.0);


    for (int i = 0; i < nParticles; i++){
        for (int j=0; j < nDimensions; j++){
            rOld(i,j) = stepLength * (distribution(gen)-0.5);
        }
    }
    rNew = rOld;

    // Create a loop over the Monte Carlo cycles here

    for (int cycle = 0; cycle < nCycles; cycle++){
        waveFuncOld = waveFunc(rOld);
        for (int i = 0; i < nParticles; i++){
            for (int j = 0; j < nDimensions; j++){
                rNew(i,j) = rOld(i,j) + stepLength*(distribution(gen)-0.5);
            }
            waveFuncNew = waveFunc(rNew);

            if (distribution(gen) <= (waveFuncNew*waveFuncNew) / (waveFuncOld*waveFuncOld)){
                for (int j = 0; j < nDimensions; j++){
                    rOld(i,j) = rNew(i,j);
                    waveFuncOld = waveFuncNew;
                }

            } else {
                for (int j = 0; j < nDimensions; j++){
                    rNew(i,j) = rOld(i,j);
                }
            }

            deltaE = localEnergy(rNew);
            energySum += deltaE;
        }
    }

    double energy = energySum/(nCycles * nParticles);
    cout << "Energy: " << energy << endl;

}



double VMC::localEnergy(const mat &r)
{
    double R1_sqrd = r(0,0)*r(0,0) + r(0,1)*r(0,1) + r(0,2)*r(0,2);
    double R2_sqrd = r(1,0)*r(1,0) + r(1,1)*r(1,1) + r(1,2)*r(1,2);
    //double R12 = sqrt(abs(R1_sqrd-R2_sqrd));
    return 0.5*omega*omega*(R1_sqrd+R2_sqrd)*(1-alpha*alpha) + 3*alpha*omega; //+ 1/R12;
}



double VMC::waveFunc(const mat &r)
{
    double R1_sqrd = r(0,0)*r(0,0) + r(0,1)*r(0,1) + r(0,2)*r(0,2);
    double R2_sqrd = r(1,0)*r(1,0) + r(1,1)*r(1,1) + r(1,2)*r(1,2);
    return exp(-alpha*omega*(R1_sqrd+R2_sqrd)*0.5);
}



