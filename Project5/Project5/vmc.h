#pragma once

#include <armadillo>

using namespace arma;

class VMC
{
public:
    VMC();
    void runMCIntegration();



private:
    double waveFunc( mat &r);
    double localEnergy( mat &r);
    double laplacian( mat &r);
    int nDimensions;
    double omega;
    double stepLength;
    int nParticles;
    double alpha;
    int nCycles;
    mat rOld;
    mat rNew;
    double accepted_states;
};
