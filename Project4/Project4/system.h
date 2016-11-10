#pragma once
#include <armadillo>
#include <string>
#include <fstream>

class System
{
public:
    System(int N_spins);
    void InitiliazeLattice(double& Energy, double& MagneticMoment);
    void MetropolisSampling(int MC_cycles, double Temperature, arma::vec &ExpectationValues, int rankProcess);
    void RunSystem(std::string filename, int MC_cycles, double initial_temp, double final_temp, double temp_step, int rankProcess, int NProcesses);
    inline int periodic(int i, int limit, int add);
    void output(std::string filename, int MC_cycles, double temperature, arma::vec expectationValues);
    //int calculateEnergy(int N_spins);
    void intermediate_output(int cycles, double temp, arma::vec expectationValues, int MC_cycles);
    //int calculate_deltaE(int ix, int iy);

private:
    arma::mat m_spinMatrix;
    int m_N_spins;
    std::ofstream m_file1;
    std::ofstream m_file2;
    double m_accepted_states;
};
