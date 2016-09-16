#pragma once
#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include "exb.cpp"
#include <algorithm>

using namespace std;


// Problem d

void exd(int argc, char* argv[]) {
    // Get the size of the matrix
    char *outfilename;
    int N;

    // Abort the program if there are to few input arguments
    if (argc <= 2)
    {
        cout << "To few aguments given" << endl;
        exit(1);
    }
    else
    {
        outfilename = argv[1];
        N = atoi(argv[2]);
        N = pow(10, N);
    }

    // -------------------------------------------------------------------------
    // Stuff

    //double h = (1./(N+1));
    double h = (1./(N+1));
    double *x = new double[N+2];
    double *b_thingy = new double[N+1];


    // Temporary arrays needed for the calculation
    double *b_hat = new double[N+1];
    double *f_hat = new double[N+1];

    // Numerical (v) and analytical (u) solution vectors
    double *u = new double[N+2];
    double *v = new double[N+2];

    //Relative error
    double *e = new double[N+2];

    // Boundary conditions
    u[0] = 0;
    v[0] = 0;
    u[N+1] = 0;
    v[N+1] = 0;

    // -------------------------------------------------------------------------

    clock_t start, finish;
    start = clock();
    // Create the x array
    for (int i = 0; i <= N+1; i++)
    {
        x[i] = i*h;
    }
    cout << x[0] << " and " << x[N+1] << endl;  // just a check!

    // We find the closed form solution
    // And we reate f*h^2, reffered to as b_thingy
    b_thingy[0] = 0;
    for (int i = 1; i <= N; i++)
    {
        b_thingy[i] = h*h*f(x[i]);
        u[i] = sol(x[i]);
    }

    // Gauss elimination and forward substitution
    b_hat[1] = 2;
    f_hat[1] = b_thingy[1];
    for (int i = 2; i <= N; i++)
    {
        b_hat[i] = (i+1.0)/i;
        f_hat[i] = b_thingy[i] + (i-1.0)*f_hat[i-1]/(i+1.0-1.0);
        //cout << f_hat[i] << " " << b_hat[i] << endl;
    }

    // Backward substitution
    v[N] = f_hat[N]/b_hat[N];
    for (int i = N-1; i >= 1; i--)
    {
        v[i] = (i)/(i+1.0)*(f_hat[i]+v[i+1]);
        //cout << v[i] << endl;
    }


    finish = clock();
    double timeused = (double) (finish - start)/(CLOCKS_PER_SEC );
    cout << setiosflags(ios::showpoint | ios::uppercase);
    cout << setprecision(10) << setw(20) << "Time used to compute our special matrix = " << timeused << " s" << endl;

    //------------------------------------------------------------------
    // Find the relative error in the data set
    double temp = 0;
    for (int i = 1; i <= N; i++)
    {
        //e[i] = abs ((v[i]-u[i])/u[i]);
        e[i] = log10 (abs ((v[i]-u[i])/u[i]));
        //cout  << v[i] << " " << u[i] << " " << e[i] << endl;
        if (abs (e[i]) > temp)
        {
            temp = e[i];
        }

    }
    cout << "The largest relative error is: " << temp << endl;

    //------------------------------------------------------------------
    // Release memory

    delete [] b_thingy;
    delete [] b_hat;
    delete [] f_hat;
    delete [] u;
    delete [] v;
    delete [] x;


}
