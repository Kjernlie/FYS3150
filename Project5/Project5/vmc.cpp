#include "vmc.h"
#include <iostream>
#include <armadillo>
#include <random>

using namespace std;



VMC::VMC() :
    nDimensions(3),
    omega(1.0),
    stepLength(1.5),
    nParticles(2),
    alpha(0.87),
    nCycles(1000000),
    accepted_states(0)
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
    energySum += localEnergy(rNew);

    // Create a loop over the Monte Carlo cycles here

    for (int cycle = 0; cycle < nCycles; cycle++){
        waveFuncOld = waveFunc(rOld);
        for (int i = 0; i < nParticles; i++){
            for (int j = 0; j < nDimensions; j++){
                rNew(i,j) = rOld(i,j) + stepLength*(distribution(gen)-0.5);
            }
            waveFuncNew = waveFunc(rNew);

            if (distribution(gen) <= (waveFuncNew*waveFuncNew) / (waveFuncOld*waveFuncOld)) {
                for (int j = 0; j < nDimensions; j++){
                    rOld(i,j) = rNew(i,j);
                }
                waveFuncOld = waveFuncNew;
                accepted_states++;

            } else {
                for (int j = 0; j < nDimensions; j++){
                    rNew(i,j) = rOld(i,j);
                }
            }

            deltaE = localEnergy(rNew);
            energySum += deltaE;
            //cout << deltaE << endl;

        }
    }

    //cout << "Energy: " << energySum << endl;
    double energy = energySum/(nCycles * nParticles);
    cout << energy << endl;
    cout << accepted_states/ (nCycles*nParticles) << endl;

}



double VMC::localEnergy(mat &r)
{
    /*
    double x1 = r(0,0);
    double y1 = r(0,1);
    double z1 = r(0,2);

    double x2 = r(1,0);
    double y2 = r(1,1);
    double z2 = r(1,2);

    double R1_sqrd = x1*x1 + y1*y1 + z1*z1;
    double R2_sqrd = x2*x2 + y2*y2 + z2*z2;
    double R12 = sqrt(  (x1-x2)*(x1-x2)  +  (y1-y2)*(y1-y2)  +  (z1-z2)*(z1-z2)  );

    */
    double R1_sqrd = r(0,0)*r(0,0) + r(0,1)*r(0,1) + r(0,2)*r(0,2);
    double R2_sqrd = r(1,0)*r(1,0) + r(1,1)*r(1,1) + r(1,2)*r(1,2);
    double R12 = sqrt( (r(0,0)-r(1,0))*(r(0,0)-r(1,0)) + (r(0,1)-r(1,1))*(r(0,1)-r(1,1)) + (r(0,2)-r(1,2))*(r(0,2)-r(1,2)) );
    return 0.5*omega*omega*(R1_sqrd + R2_sqrd)*(1 - alpha*alpha) + 3*alpha*omega + 1./R12;


    /*
    double R1_sqrd = r(0,0)*r(0,0) + r(0,1)*r(0,1)+ r(0,2)*r(0,2);
    double R2_sqrd = r(1,0)*r(1,0) + r(1,1)*r(1,1)+ r(1,2)*r(1,2);
    double L       = laplacian(r);
    double R12 = sqrt( (r(0,0)-r(1,0))*(r(0,0)-r(1,0)) + (r(0,1)-r(1,1))*(r(0,1)-r(1,1)) + (r(0,2)-r(1,2))*(r(0,2)-r(1,2)) );
    return 0.5*omega*omega*(R1_sqrd + R2_sqrd) + L + 1/R12;
    */
}



double VMC::waveFunc( mat &r)
{
    double R1_sqrd = r(0,0)*r(0,0) + r(0,1)*r(0,1) + r(0,2)*r(0,2);
    double R2_sqrd = r(1,0)*r(1,0) + r(1,1)*r(1,1) + r(1,2)*r(1,2);
    return exp(-alpha*omega*(R1_sqrd+R2_sqrd)*0.5);
}

//double VMC::laplacian(mat &r) {
//    double L = 0;
//    double wf = waveFunc(r);
//    double dx = 1e-5;

//    for (int i = 0; i < nParticles; i++) {
//        for (int j = 0; j < nDimensions; j++) {
//            r(i,j) += dx;
//            double wfp = waveFunc(r);
//            r(i,j) -= 2*dx;
//            double wfm = waveFunc(r);
//            r(i,j) += dx;
//            L += wfp+wfm-2*wf;
//        }
//    }
//    L *= -0.5 / (dx*dx * wf);
//    return L;
//}


