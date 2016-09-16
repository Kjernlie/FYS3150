#include <iostream>
#include <armadillo>
#include <fstream>
using namespace arma;
using namespace std;


// Problem e)


void exe(int argc, char* argv[])
{
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

    // Start time
    clock_t start, finish;
    start = clock();

    // Create h, define some vectors and matrix A
    double h = 1./(N+1);
    vec x = linspace<vec>(0, 1, N+2);
    vec f = vec(N);
    mat A = mat(N, N, fill::zeros);
    vec b = vec(N);

    // Fill the right hand side vector
    for (int i = 0; i <= N-1; i++)
    {
        f(i) = 100*exp(-10*x(i+1))*h*h;
    }

    // Fill the tridiagonal matrix
    A.diag(1).fill(-1);
    A.diag().fill(2);
    A.diag(-1).fill(-1);

    // Solve with LU decomposition
    mat L, U, P;
    lu(L,U,P,A);
    mat F = P * f;

    vec y = solve(L,F);
    vec z = solve(U,y);

    // Stop time
    finish = clock();
    double timeused = (double) (finish - start)/(CLOCKS_PER_SEC );
    cout << "Time used with LU decomposition = " << timeused << " s" << endl;

    // The build-in function solve will also solve the system with LU decomposition
    vec v = solve(A,f);
    cout << v << endl;



    // ------------------------------------------------------------------
    // Write to file

    ofstream myfile;
    myfile.open(outfilename);
    //myfile << setiosflags(ios::showpoint | ios::uppercase);
    for (int i = 0; i <= N-1; i++)
    {
        myfile << v[i] << endl;
    }
    myfile.close();

}

