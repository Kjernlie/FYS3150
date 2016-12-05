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



VMC::VMC(int T) :
    nDimensions(3),
    omega(1.0),
    stepLength(1.0),
    nParticles(2),
    alpha(0.9),
    nCycles(2000000),
    accepted_states(0),
    h(0.00001),
    h2(1./double(h*h)),
    beta(0.05),
    perturbation(1),
    m_T(T),
    m_r12(0)
{
}


void VMC::runMCIntegration()
{

    ostringstream oss;
    oss << "data_w" << omega << "_N_" << nCycles << ".dat";
    string filename = oss.str();

    m_file.open(filename);


    for (alpha = 0.6; alpha <= 1.2; alpha += 0.05)
    {

        //optimalAlpha();
        //optimalAlphaBeta();

        accepted_states = 0;

        rOld = zeros<mat>(nParticles, nDimensions);
        rNew = zeros<mat>(nParticles, nDimensions);
        double energySum = 0;
        double energySum2 = 0;
        double meanDistance = 0;
        MCIntegration(energySum, energySum2, meanDistance);

        double energy = energySum/(nCycles * nParticles);
        double energy2 = energySum2/(nCycles * nParticles);
        double variance = energy2 - energy*energy;
        cout << "Energy: " << energy << endl;
        cout << "Variance: " << variance << endl;
        cout << "Acceptance rate: " << accepted_states/ (nCycles*nParticles) << endl;
        cout << "Step length: " << stepLength << endl;
        cout << "Mean distance: " << meanDistance/((double)nCycles/1000) << endl;

        writeToFile(energy, variance);
    }

}



void VMC::MCIntegration(double& energySum, double& energySum2, double &meanDistance)
{
    // Initialize the seed and calle the Mersenne algorithm
    random_device rd;
    mt19937_64 gen(rd());

    // Setup the uniform distribution for x in [0,1]
    uniform_real_distribution<double> distribution(0.0,1.0);


    double waveFuncOld = 0;
    double waveFuncNew = 0;
    double deltaE = 0;


    for (int i = 0; i < nParticles; i++){
        for (int j=0; j < nDimensions; j++){
            rOld(i,j) = stepLength * (distribution(gen)-0.5);
        }
    }
    rNew = rOld;
    deltaE = localEnergy(rNew);
    energySum += deltaE;
    energySum2 += deltaE*deltaE;


    // Loop over MC cycles
    for (int cycle = 0; cycle < nCycles; cycle++){
        waveFuncOld = (m_T == 1) ? waveFunc(rOld) : waveFunc2(rOld);
        for (int i = 0; i < nParticles; i++){
            for (int j = 0; j < nDimensions; j++){
                rNew(i,j) = rOld(i,j) + stepLength*(distribution(gen)-0.5);
            }
            waveFuncNew = (m_T == 1) ? waveFunc(rNew) : waveFunc2(rNew);

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


        }



        // Adjust the step length towards an optimal value
        if (cycle % 1000 == 0 && cycle < 10000){
            acceptanceTest(cycle);
        }
        if (cycle % 1000 == 0){
            //meanDistance += rDiff(rNew);
            meanDistance += m_r12;
        }

    }
}




double VMC::localEnergy(const mat &r)
{
    mat rPlus = zeros<mat>(nParticles, nDimensions);
    mat rMinus = zeros<mat>(nParticles, nDimensions);
    rPlus = rMinus = r;

    double waveFuncPlus = 0;
    double waveFuncMinus = 0;
    double waveFuncCurrent = (m_T == 1) ? waveFunc(r) : waveFunc2(r);

    double kineticEnergy = 0;

    for (int i = 0; i < nParticles; i++){
        for (int j = 0; j < nDimensions; j++){
            rPlus(i,j) += h;
            rMinus(i,j) -= h;
            waveFuncMinus = (m_T == 1) ? waveFunc(rMinus) : waveFunc2(rMinus);
            waveFuncPlus = (m_T == 1) ? waveFunc(rPlus) : waveFunc2(rPlus);
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
                m_r12 = sqrt(r12);
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
            argument2 += sqrt(r12);
        }
    }

    return exp(-0.5 * argument * alpha * omega) * exp(argument2/ ( 2*( 1 + beta*argument2) ) );
}


double VMC::rDiff(const mat &r)
{

    double r12 = 0;
    double argument2 = 0;
    for (int i = 0; i < nParticles; i++){
        for (int j = i+1; j < nParticles; j++){
            r12 = 0;
            for (int k = 0; k < nDimensions; k++){
                r12 += (r(i,k) - r(j,k)) * (r(i,k) - r(j,k));
            }
            argument2 += sqrt(r12);
        }
    }

    return argument2;
}


void VMC::acceptanceTest(int cycles){
    double accepted_prob = accepted_states/(cycles*nParticles);
    if (accepted_prob > 0.6){
        stepLength += 0.1;
    } else if (accepted_prob < 0.4){
        stepLength -= 0.1;
    }
}


void VMC::optimalAlpha()
{
    double alphaMax = 0.72;
    double alphaMin = 0.3;
    double stepSize = 0.02;

    int N = (int) ((alphaMax - alphaMin) / stepSize);

    vec energies = zeros<vec>(N);
    vec alphas = zeros<vec>(N);

    int counter = 0;
    for (alpha = alphaMin; alpha <= alphaMax; alpha += stepSize){

        accepted_states = 0;
        rOld = zeros<mat>(nParticles, nDimensions);
        rNew = zeros<mat>(nParticles, nDimensions);
        double energySum = 0;
        double energySum2 = 0;
        double meanDistance = 0;
        MCIntegration(energySum, energySum2, meanDistance);
        energies(counter) = energySum;
        alphas(counter) = alpha;
        counter++;

    }

    double min = energies(0);
    int minIndex = 0;
    for (int i = 0; i < N; i++){
        if (energies(i) < min){
            min = energies(i);
            minIndex = i;
        }
    }

    alpha = alphas(minIndex);

    cout << "Optimal alpha: " << alpha << endl;

}


void VMC::optimalAlphaBeta()
{

    double alphaMax = 1.1;
    double alphaMin = 0.7;
    double betaMax = 0.7;
    double betaMin = 0.0;
    double stepSize = 0.05;


    //int N = (int) ((alphaMax - alphaMin)/stepSize) * ((betaMax - betaMin)/stepSize);

    int N = ((alphaMax - alphaMin) / stepSize) * ((betaMax - betaMin) / stepSize);
    cout << N << endl;

    vec energies = zeros<vec>(N);
    vec alphas = zeros<vec>(N);
    vec betas = zeros<vec>(N);


    int counter = 0;
    for (alpha = alphaMin; alpha <= alphaMax; alpha += stepSize){
        for (beta = betaMin; beta <= betaMax; beta += stepSize){
            accepted_states = 0;
            rOld = zeros<mat>(nParticles, nDimensions);
            rNew = zeros<mat>(nParticles, nDimensions);
            double energySum = 0;
            double energySum2 = 0;
            double meanDistance = 0;
            MCIntegration(energySum, energySum2, meanDistance);
            energies(counter) = energySum;
            alphas(counter) = alpha;
            betas(counter) = beta;
            counter++;
        }
    }

    double min = energies(0);
    int minIndex = 0;
    for (int i = 0; i < N; i++){
        if (energies(i) < min){
            min = energies(i);
            minIndex = i;
        }
    }

    alpha = alphas(minIndex);
    beta = betas(minIndex);

    cout << "Optimal alpha: " << alpha << endl;
    cout << "Optimal beta: " << beta << endl;


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
