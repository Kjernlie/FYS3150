#pragma once

#include <armadillo>

using namespace arma;

class VMC
{
public:
    VMC();
    void runMCIntegration();



private:
    double waveFunc(const mat &r);
    double localEnergy(const mat &r);
    int nDimensions;
    int omega;
    double stepLength;
    int nParticles;
    double alpha;
    int nCycles;
    mat rOld;
    mat rNew;
};
