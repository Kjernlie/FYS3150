#include <iostream>
#include "vmc.h"


using namespace std;

int main()
{
    VMC solver;
    solver.runMCIntegration();

    cout << "So, this works?" << endl;
    return 0;
}

