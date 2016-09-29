#include <iostream>
#include <armadillo>
#include "jacobi.h"
#include "exb.h"
#include <fstream>

using namespace std;
using namespace arma;

void exb()
{

//    mat BenchMrk = mat(3,3,fill::zeros);
//    BenchMrk.diag().fill(1);
//    BenchMrk(0,1) = 2;
//    BenchMrk(1,0) = 2;

//    cout << BenchMrk << endl;
//    mat R = mat(3,3,fill::zeros);
//    R.diag().fill(1);

//    int counter = 0;
//    vec max;
//    while (testing(BenchMrk) == 0)
//    {
//        counter++;
//        max = max_element(BenchMrk);
//        rotate(BenchMrk,R,max(1),max(2));

//    }

//    vec eig_vals = BenchMrk.diag();
//    eig_vals = sort(eig_vals);

//    cout << "first eigenvalue: " << eig_vals(0) << endl;
//    cout << "second eigenvalue: " << eig_vals(1) << endl;
//    cout << "third eigenvalue: " << eig_vals(2) << endl;
//    cout << "number of similarity transformations: " << counter << endl;


    // Initialize parameters
    double rho_min = 0;
    double rho_max = 10;
    int N = 100;

    // Initialize constants
    double h = rho_max/ (N+1);
    double diag_const = 2.0/(h*h);
    double nondiag_const = -1.0/(h*h);

    // Calculate array of potential values

    vec v(N);
    vec rho(N);
    rho(0) = rho_min;
    v(0) = rho(0)*rho(0);
    for (int i = 0; i<N; i++)
    {
        rho(i) = rho_min + (i+1)*h;
        v(i) = rho(i)*rho(i);
    }


    // Creating the threedimensional matrix
    mat A = mat(N,N,fill::zeros);
    A.diag(1).fill(nondiag_const);
    A.diag(-1).fill(nondiag_const);
    A.diag().fill(diag_const);


    for (int i = 0; i < N; i++)
    {
        //double rho = (i+1)*h;
        //A(i,i) = A(i,i) + rho*rho;
        A(i,i) = A(i,i) + v(i);
    }

    // eigenvector
    mat R = mat(N,N,fill::zeros);
    R.diag().fill(1);

    int counter = 0;

    vec max;
    while (testing(A) == 0)
    {
        counter++;
        max = max_element(A);
        rotate(A,R,max(1),max(2));

    }

    //cout << R << endl;
    vec eig_vals = A.diag();
    uvec indices = sort_index(eig_vals);
    eig_vals = sort(eig_vals);


    cout << eig_vals(0) << endl;
    cout << eig_vals(1) << endl;
    cout << eig_vals(2) << endl;
    cout << "number of similarity transformations: " << counter << endl;




    // ------------------------------------------------------------------------------------------------

   // Write to file
    ofstream myfile;
    myfile.open("data.txt");

    for (int i = 0; i < N; i++)
    {
        myfile << rho(i) << " " << R(i,indices(0)) << " " << R(i,indices(1)) << " " << R(i,indices(2)) << endl;
    }
    myfile.close();

}
