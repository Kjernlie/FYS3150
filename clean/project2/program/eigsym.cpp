#include <iostream>
#include <armadillo>
#include <eigsym.h>
#include <iomanip>

using namespace std;
using namespace arma;

void eigsym()
{
    // Initialize parameters
    int rho_min = 0;
    int rho_max = 5;
    int N = 400;

    // Initialize constants
    double h = double(rho_max)/(N);
    double diag_const = 2.0/(h*h);
    double nondiag_const = -1.0/(h*h);

    mat A = mat(N,N,fill::zeros);
    A.diag(1).fill(nondiag_const);
    A.diag(-1).fill(nondiag_const);
    A.diag().fill(diag_const);


    for (int i = 0; i < N; i++)
    {
        double rho = (i+1)*h;
        A(i,i) = A(i,i) + rho*rho;
    }

    clock_t start, finish;
    start = clock();

    // using the built-in function
    vec eigval;
    mat eigvec;
    eig_sym(eigval,eigvec,A);
    cout << eigval(0) << endl;
    cout << eigval(1) << endl;
    cout << eigval(2) << endl;


    // -----------------------------------------------------------------------------------
    // Time

    finish = clock();
    double timeused = (double) (finish - start)/(CLOCKS_PER_SEC );
    cout << "Time used for a threediagonal matrix = " << timeused << " s" << endl;

    // -------------------------------------------------------------------------------------
//    for (int i = 0; i < N; i++)
//    {
//        cout << eigvec(0,i)
//    }
//    cout << eigvec.n_cols << endl;
}
