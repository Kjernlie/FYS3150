#include <iostream>
#include "vmc.h"



using namespace std;

int main()
{
    clock_t start, finish;
    start = clock();

    int T = 1;
    VMC solver(T);
    solver.runMCIntegration();

    finish = clock();
    double timeused = (double) (finish - start)/(CLOCKS_PER_SEC );
    cout << "Time used is  = " << timeused << " s" << endl;

    cout << "So, this works?" << endl;
    return 0;
}

