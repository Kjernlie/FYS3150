#include <iostream>
#include "vmc.h"
double ran2(long *);

using namespace std;

int main()
{
    VMC solver;
    solver.runMCIntegration();
    return 0;
}
