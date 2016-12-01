#ifndef VMC_H
#define VMC_H
#include <armadillo>

using namespace arma;

class VMC
{
public:
    VMC();
    runMCIntegration;

private:
    double waveFunc(const mat &r);
    double localEnergy(const mat &r);
    int nDimensions;
    int omega;
    double stepLength;
    int nParticles;
    long idum;
    double alpha;
    int nCycles;
    mat rOld;
    mat rNew;
};

#endif // VMC_H
