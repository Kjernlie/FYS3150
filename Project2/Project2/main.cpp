#include <iostream>
#include <jacobi.h>
#include <armadillo>

using namespace std;
using namespace arma;

// Project 2

int N = 2;
mat A = mat(N, N, fill::ones);


int main()
{

    vec test = vec(3,fill::zeros);
    vec test2 = vec(2,fill::zeros);
    mat test3 = mat(N,N,fill::zeros);

    test = max_element(A);
    test2 = trigonometry(A,test(1),test(2));
    test3 = make_B(A,test(1),test(2),test2(0),test(1));
    cout << test3 << endl;

    //cout << max_element(A) << endl;
    return 0;
}


