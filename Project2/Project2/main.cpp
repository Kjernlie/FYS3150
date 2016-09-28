#include <iostream>
#include <jacobi.h>
#include <armadillo>

using namespace std;
using namespace arma;

// Project 2



int main()
{

    // Initialize parameters
    int rho_min = 0;
    int rho_max = 10;
    int N = 20;

    // Initialize constants
    double h = double(rho_max)/N;
    double diag_const = 2.0/(h*h);
    double nondiag_const = -1.0/(h*h);

    // Calculate array of potential values
    vec v(N);
    vec rho(N);
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
        A(i,i) = A(i,i) + v(i);
    }
    //cout << A << endl;


//    // using the built-in function
//    vec eigval;
//    mat eigvec;
//    eig_sym(eigval,eigvec,A);
//    cout << eigval(0) << endl;
//    cout << eigval(1) << endl;
//    cout << eigval(2) << endl;

    mat B;
    vec max;
    vec trig;
    //int tol = 1e-8;
    int counter = 0;


    //default tolerance value is set to 1e-8 in header file
    while (testing(A) == 0)
    {
        counter++;
        max = max_element(A);
        trig = trigonometry(A,max(1),max(2));
        B = make_B(A,max(1),max(2),trig(0),trig(1));
        A = B;
    }


    cout << counter << endl;
    cout << A(0,0) << endl;
    cout << A(1,1) << endl;
    cout << A(2,2) << endl;


    return 0;
}


