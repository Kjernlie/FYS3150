#pragma once
#include <armadillo>
#include <fstream>

using namespace arma;

class VMC
{
public:
    VMC();
    void runMCIntegration();
    void writeToFile(double energy, double variance);



private:
    double waveFunc(const mat &r);
    double localEnergy(const mat &r);
    void acceptanceTest(int cycles);
    int nDimensions;
    double omega;
    double stepLength;
    int nParticles;
    double alpha;
    int nCycles;
    mat rOld;
    mat rNew;
    double accepted_states;
    double h;
    double h2;
    int perturbation;
    std::ofstream m_file;
};
