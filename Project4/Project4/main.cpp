#include <iostream>
#include "system.h"

using namespace std;

int main()
{
    int N_spins = 2;
    int MC_cycles = 1e5;
    double initial_temp = 0.0;
    double final_temp = 1.0;
    double temp_step = 0.1;

    System test;

    test.RunSystem(MC_cycles, N_spins, initial_temp, final_temp, temp_step);



    cout << "Hello World!" << endl;
    return 0;
}
