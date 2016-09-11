#include <iostream>
#include <armadillo>
using namespace arma;
using namespace std;


// Problem 2 e)


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


    double h = 1./(N+1);
    vec x = linspace<vec>(0, 1, N+2);
    vec f = vec(N);
    mat A = mat(N, N, fill::zeros);
    vec b = vec(N);

    for (int i = 0; i <= N-1; i++)
    {
        f(i) = 100*exp(-10*x(i+1))*h*h;
    }

    A.diag(1).fill(-1);
    A.diag().fill(2);
    A.diag(-1).fill(-1);


    mat L, U, P;
    lu(L,U,P,A);
    cout << P << endl;
    mat F = P * f;
    cout << F << endl;

    vec y = solve(L,F);
    vec z = solve(U,y);

    vec v = solve(A,f);
    cout << v-z << endl;

}
