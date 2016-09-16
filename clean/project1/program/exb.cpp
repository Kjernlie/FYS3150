#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include "exb.h"


using namespace std;

// Problem b)


// create source function
double f(double x)
{
    return 100*exp(-10*x);
}
// closed form solution
double sol(double x)
{
    return 1 - (1 - exp(-10))*x - exp(-10*x);
}

// ----------------------------------------------------------------------------

void exb(int argc, char* argv[]) {
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
    }

    // -------------------------------------------------------------------------
    // Stuff

    clock_t start, finish;
    start = clock();

    double h = 1./(N+1);
    double *x = new double[N+2];
    double *b_thingy = new double[N+1];

    // Temporary arrays needed for the calculation
    double *b_hat = new double[N+1];
    double *f_hat = new double[N+1];


    // Numerical (v) and analytical (u) solution vectors
    double *u = new double[N+2];
    double *v = new double[N+2];

    // Boundary conditions
    u[0] = 0;
    v[0] = 0;
    u[N+1] = 0;
    v[N+1] = 0;

    // -------------------------------------------------------------------------


    // Create the arrays in the three diagonal matrix
    // using dynamic memory allocation
    int *a = new int[N+1];
    int *b = new int[N+1];
    int *c = new int[N+1];



    // Create the x array
    for (int i = 0; i <= N+1; i++)
    {
        x[i] = i*h;
    }
    cout << x[0] << " and " << x[N+1] << endl;  // just a check!

    // Create f*h^2, reffered to as b_thingy
    // Also, we define the values on the diagonals
    // Furthermore, we find the closed form solution
    b_thingy[0] = 0;
    for (int i = 1; i <= N; i++)
    {
        b_thingy[i] = h*h*f(x[i]);

        // values on the diagonals
        a[i] = -1.0;
        b[i] = 2.0;
        c[i] = -1.0;

        // closed form solution
        u[i] = sol(x[i]);
    }
    a[1] = 0;
    c[N] = 0;

    // Forward substitution
    b_hat[1] = b[1];
    f_hat[1] = b_thingy[1];
    for (int i = 2; i <= N; i++)
    {
        b_hat[i] = b[i] - (a[i]*c[i-1])/b_hat[i-1];
        f_hat[i] = b_thingy[i] - (a[i]*f_hat[i-1])/b_hat[i-1];
    }

    // Backward substitution
    v[N] = f_hat[N]/b_hat[N];
    for (int i = N-1; i >= 1; i--)
    {
        v[i] = (f_hat[i]-c[i]*v[i+1])/b_hat[i];
        cout << v[i] << endl;
    }


    // --------------------------------------------------------------------
    // Stop time, and print to screen

    finish = clock();
    double timeused = (double) (finish - start)/(CLOCKS_PER_SEC );
    cout << "Time used for a threediagonal matrix = " << timeused << " s" << endl;


    // ----------------------------------------------------------------
    // Write to file

    ofstream myfile;
    myfile.open(outfilename);
    //myfile.open("balle.dat");
    //myfile << setiosflags(ios::showpoint | ios::uppercase);
    for (int i = 0; i <= N+1; i++)
    {
        myfile << x[i] << " " << u[i] << " " << v[i] << endl;
    }
    myfile.close();


    // --------------------------------------------------------------------
    // release memory


    delete [] a;
    delete [] b;
    delete [] c;
    delete [] b_thingy;
    delete [] b_hat;
    delete [] f_hat;
    delete [] u;
    delete [] v;
    delete [] x;

}
