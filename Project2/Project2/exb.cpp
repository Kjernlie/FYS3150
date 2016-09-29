#include <iostream>
#include <armadillo>
#include "jacobi.h"
#include "exb.h"

using namespace std;
using namespace arma;

void exb()
{


    // Initialize parameters
    double rho_min = 0;
    double rho_max = 6;
    int N = 100;

    // Initialize constants
    double h = rho_max/ (N+1);
    double diag_const = 2.0/(h*h);
    double nondiag_const = -1.0/(h*h);

    // Calculate array of potential values
    /*
    vec v = zeros<vec>(N);
    vec rho(N);
    rho(0) = rho_min;
    v(0) = rho(0)*rho(0);
    for (int i = 0; i<N; i++)
    {
        rho(i) = rho_min + (i+1)*h;
        v(i) = rho(i)*rho(i);
    }
    */

    // Creating the threedimensional matrix
    mat A = mat(N,N,fill::zeros);
    A.diag(1).fill(nondiag_const);
    A.diag(-1).fill(nondiag_const);
    A.diag().fill(diag_const);


    for (int i = 0; i < N; i++)
    {
        double rho = (i+1)*h;
        A(i,i) = A(i,i) + rho*rho;
    }

    int counter = 0;

    vec max;
    while (testing(A) == 0)
    {
        counter++;
        max = max_element(A);
        rotate(A,max(1),max(2));

    }

    cout << counter << endl;
    vec eig_vals = A.diag();
    eig_vals = sort(eig_vals);

    cout << eig_vals(0) << endl;
    cout << eig_vals(1) << endl;
    cout << eig_vals(2) << endl;
    cout << "number of similarity transformations: " << counter << endl;

}
