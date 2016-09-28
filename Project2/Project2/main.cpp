#include <iostream>
#include <jacobi.h>
#include <armadillo>

using namespace std;
using namespace arma;

// Project 2



int main()
{

//    int rho_min = 0.0;
//    int rho_max = 10.0;
//    int lOrbital = 0;
//    int dim = 400;


    int N = 2;
    mat A = mat(N,N,fill::ones);
    A(0,0) = 2;
    A(1,1) = 2;

    cout << A << endl;

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
        cout << counter << endl;
        cout << A << endl;
    }

    return 0;
}


