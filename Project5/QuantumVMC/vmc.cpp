#include "vmc.h"
#include <iostream>
#include <armadillo>
#include <random>
#include <string>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>

using namespace std;



VMC::VMC() :
    nDimensions(3),
    omega(1.0),
    stepLength(1.0),
    nParticles(2),
    alpha(1.0),
    nCycles(100000),
    accepted_states(0),
    h(0.00001),
    h2(1./double(h*h)),
    beta(0.6),
    perturbation(1)
{
}


void VMC::runMCIntegration()
{
    // Initialize the seed and calle the Mersenne algorithm
    random_device rd;
    mt19937_64 gen(rd());

    // Setup the uniform distribution for x in [0,1]
    uniform_real_distribution<double> distribution(0.0,1.0);

    ostringstream oss;
    oss << "testing_w" << omega << "_N_" << nCycles << ".dat";
    string filename = oss.str();

    m_file.open(filename);


    for (alpha = 0.9; alpha <= 1.1; alpha += 0.1)
    {


        accepted_states = 0;
        rOld = zeros<mat>(nParticles, nDimensions);
        rNew = zeros<mat>(nParticles, nDimensions);
        double waveFuncOld = 0;
        double waveFuncNew = 0;
        double energySum = 0;
        double energySum2 = 0;
        double deltaE;


        for (int i = 0; i < nParticles; i++){
            for (int j=0; j < nDimensions; j++){
                rOld(i,j) = stepLength * (distribution(gen)-0.5);
            }
        }
        rNew = rOld;
        //energySum += localEnergy(rNew);

        // Create a loop over the Monte Carlo cycles here

        for (int cycle = 0; cycle < nCycles; cycle++){
            waveFuncOld = waveFunc2(rOld);
            for (int i = 0; i < nParticles; i++){
                for (int j = 0; j < nDimensions; j++){
                    rNew(i,j) = rOld(i,j) + stepLength*(distribution(gen)-0.5);
                }
                waveFuncNew = waveFunc2(rNew);

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
                energySum2 += deltaE*deltaE;
                //cout << deltaE << endl;

            }

            // Find the optimal value for the step length
            if (cycle % 1000 == 0 && cycle < 10000){
                acceptanceTest(cycle);
            }

        }

        double energy = energySum/(nCycles * nParticles);
        double energy2 = energySum2/(nCycles * nParticles);
        double variance = energy2 - energy*energy;
        cout << "Energy: " << energy << endl;
        cout << "Variance: " << variance << endl;
        cout << "Acceptance rate: " << accepted_states/ (nCycles*nParticles) << endl;

        writeToFile(energy, variance);
    }

}



double VMC::localEnergy(const mat &r)
{
    mat rPlus = zeros<mat>(nParticles, nDimensions);
    mat rMinus = zeros<mat>(nParticles, nDimensions);
    rPlus = rMinus = r;

    double waveFuncPlus = 0;
    double waveFuncMinus = 0;
    double waveFuncCurrent = waveFunc2(r);

    double kineticEnergy = 0;

    for (int i = 0; i < nParticles; i++){
        for (int j = 0; j < nDimensions; j++){
            rPlus(i,j) += h;
            rMinus(i,j) -= h;
            waveFuncMinus = waveFunc2(rMinus);
            waveFuncPlus = waveFunc2(rPlus);
            kineticEnergy -= (waveFuncMinus + waveFuncPlus - 2 * waveFuncCurrent);
            rPlus(i,j) = r(i,j);
            rMinus(i,j) = r(i,j);

        }
    }

    kineticEnergy *= 0.5 * h2 / waveFuncCurrent;

    double potentialEnergy = 0;
    double r_oneParticle = 0;


    for (int i = 0; i < nParticles; i++){
        r_oneParticle =0;
        for (int j = 0; j < nDimensions; j++){
            r_oneParticle += r(i,j) * r(i,j);
        }

        potentialEnergy += r_oneParticle;
    }

    potentialEnergy *= 0.5 * omega * omega;

    if (perturbation == 1){
        double r12 = 0;
        for (int i = 0; i < nParticles; i++){
            r_oneParticle = 0;
            for (int j = i+1; j < nParticles; j++){
                r12 = 0;
                for (int k = 0; k < nDimensions; k++){
                    r12 += (r(i,k) - r(j,k)) * (r(i,k) - r(j,k));
                }
                potentialEnergy += 1/sqrt(r12);
            }
        }
    }

    return kineticEnergy + potentialEnergy;


}




double VMC::waveFunc(const mat &r){
    double argument = 0;
    for (int i = 0; i < nParticles; i++){
        double r_oneParticle = 0;
        for (int j = 0; j < nDimensions; j++){
            r_oneParticle += r(i,j) * r(i,j);
        }

        argument += r_oneParticle;
    }

    return exp(-0.5 * argument * alpha * omega);
}


double VMC::waveFunc2(const mat &r)
{
    double argument = 0;
    for (int i = 0; i < nParticles; i++){
        double r_oneParticle = 0;
        for (int j = 0; j < nDimensions; j++){
            r_oneParticle += r(i,j) * r(i,j);
        }

        argument += r_oneParticle;
    }

    double r12 = 0;
    double argument2 = 0;
    for (int i = 0; i < nParticles; i++){
        for (int j = i+1; j < nParticles; j++){
            r12 = 0;
            for (int k = 0; k < nDimensions; k++){
                r12 += (r(i,k) - r(j,k)) * (r(i,k) - r(j,k));
            }
            argument2 += r12;
        }
    }

    return exp(-0.5 * argument * alpha * omega)*exp(argument2/ ( 2*( 1 + beta*argument2) ) );
}


void VMC::acceptanceTest(int cycles){
    double accepted_prob = accepted_states/(cycles*nParticles);
    if (accepted_prob > 0.6){
        stepLength += 0.1;
    } else if (accepted_prob < 0.4){
        stepLength -= 0.1;
    }
}


void VMC::writeToFile(double energy, double variance){


// Why isn't the file good?

//    if(!m_file.good()) {
//        m_file.open(filename.c_str(), ofstream::out);
//        if(!m_file.good()) {
//            cout << "Error opening file " << filename << ". Aborting!" << endl;
//            terminate();
//        }
//    }


    m_file << setiosflags(ios::showpoint | ios::uppercase);

    m_file << setw(15) << setprecision(8) << alpha;
    m_file << setw(15) << setprecision(8) << energy;
    m_file << setw(15) << setprecision(8) << variance;
    m_file << "\n";


}
