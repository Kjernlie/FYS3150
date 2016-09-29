#include <iostream>
#include <armadillo>
#include <eigsym.h>

using namespace std;
using namespace arma;

void eigsym()
{
    // Initialize parameters
    int rho_min = 0;
    int rho_max = 10;
    int N = 400;

    // Initialize constants
    double h = double(rho_max)/(N);
    double diag_const = 2.0/(h*h);
    double nondiag_const = -1.0/(h*h);

    // Calculate array of potential values
//    vec v(N);
//    vec rho(N);
//    rho(0) = rho_min;
//    v(0) = rho(0)*rho(0);
//    for (int i = 0; i<N; i++)
//    {
//        rho(i) = rho_min + (i+1)*h;
//        v(i) = rho(i)*rho(i);
//    }

//    // Creating the threedimensional matrix
//    mat A = mat(N,N,fill::zeros);
//    A.diag(1).fill(nondiag_const);
//    A.diag(-1).fill(nondiag_const);
//    A.diag().fill(diag_const);


//    for (int i = 0; i < N; i++)
//    {
//        A(i,i) = A(i,i) + v(i);
//    }

    mat A = mat(N,N,fill::zeros);
    A.diag(1).fill(nondiag_const);
    A.diag(-1).fill(nondiag_const);
    A.diag().fill(diag_const);


    for (int i = 0; i < N; i++)
    {
        double rho = (i+1)*h;
        A(i,i) = A(i,i) + rho*rho;
    }


    // using the built-in function
    vec eigval;
    mat eigvec;
    eig_sym(eigval,eigvec,A);
    cout << eigval(0) << endl;
    cout << eigval(1) << endl;
    cout << eigval(2) << endl;
//    for (int i = 0; i < N; i++)
//    {
//        cout << eigvec(0,i)
//    }
//    cout << eigvec.n_cols << endl;
}
