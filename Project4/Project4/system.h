#pragma once
#include <armadillo>
#include <string>
#include <fstream>

class System
{
public:
    System();
    void InitiliazeLattice(int N_spins, arma::mat &SpinMatrix, double& Energy, double& MagneticMoment);
    void MetropolisSampling(int N_spins, int MC_cycles, double Temperature, arma::vec &ExpectationValues);
    void RunSystem(std::string filename, int MC_cycles, int N_spins, double initial_temp, double final_temp, double temp_step, int rankProcess, int NProcesses);
    inline int periodic(int i, int limit, int add);
    void output(std::string filename, int N_spins, int MC_cycles, double temperature, arma::vec ExpectationValues);
    int calculateEnergy(arma::mat spinMatrix, int N_spins);

private:
    std::ofstream m_file;
};
