#pragma once
#include <armadillo>
#include <fstream>

using namespace arma;

class VMC
{
public:
    VMC(int T);
    void runMCIntegration();
    void writeToFile(double energy, double variance);



private:
    void MCIntegration(double &energySum, double &energySum2, double &meanDistance);
    double waveFunc(const mat &r);
    double waveFunc2(const mat &r);
    double localEnergy(const mat &r);
    void acceptanceTest(int cycles);
    void optimalAlphaBeta();
    void optimalAlpha();
    double rDiff(const mat &r);
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
    double beta;
    double m_r12;
    int m_T;
};
