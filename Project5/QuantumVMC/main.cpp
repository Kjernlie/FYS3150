#include <iostream>
#include "vmc.h"
#include <armadillo>



using namespace std;

int main()
{
    clock_t start, finish;
    start = clock();

    int T1 = 1;
    int T2 = 2;
    double omegas [3] = {0.01, 0.5, 1.0};
    double stepLengths [3] = {10, 1.4, 1.0};

    VMC solver(T1);

    for (int i = 0; i < 3; i++){
        solver.omega = omegas[i];
        cout << "Omega: " << omegas[i] << endl;
        solver.stepLength = stepLengths[i];
        solver.runMCIntegration();
        cout << "\n" << endl;
    }


    finish = clock();
    double timeused = (double) (finish - start)/(CLOCKS_PER_SEC );
    cout << "Time used is  = " << timeused << " s" << endl;

    cout << "So, this works?" << endl;
    return 0;
}

