#include "jacobi.h"
#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

// What to do!

// 1. choose a tolerance....this is maybe better to do in main
// 2. while test... compare norms...also maybe main, unless....
// 3. choose the largest matrix element akl...create a function here
// 4. compute tau, sin, cos and tan...create a function here
// 5. compute the similarity transformation....create a function here
// 6. compute the new norm of the off-diagonal elements.... create function here
// 7. and cont. till norm(B) < is satisfied....this is the second point


// ------------------------------------------------------------------------------------------------------------------
// Find the largest non-diagonal element in A!

vec max_element(mat A)
{
    int N = A.n_cols;
    vec max = vec(3,fill::zeros); // max(0) is the max value, max(1) and max(2) will represent the indices, k and l, respectively

    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            if (i != j)
            {
                double aij = abs(A(i,j));
                if (aij > max(0))
                {
                    max(0) = aij, max(1) = i, max(2) = j;
                }

            }
        }
    }

    return max;
}

// ----------------------------------------------------------------------------------------------------------------------
// Find the trigonometric elements


vec trigonometry(mat A, int k, int l)
{
    // Do I have to put in a case for akl = 0?

    vec trig = vec(2,fill::zeros);
    double tau, t;

    tau = (A(l,l) - A(k,k))/(2*A(k,l));

    if (tau >= 0)
    {
        t = 1.0/(tau + sqrt(1.0+tau*tau));
    }
    else
    {
        t = -1.0/(-tau + sqrt(1.0+tau*tau));
    }

    trig(0) = 1.0/sqrt(1+t*t);          // cosine
    trig(1) = t*trig(0);                // sine

    return trig;
}



// ---------------------------------------------------------------------------------------------------------------------------
// create the new matrix B

mat make_B(mat A, int k, int l, double c, double s)
{
    int N = A.n_cols;               // find the mesh size
    mat B = mat(N,N,fill::zeros);   // create matrix B filled with zeros

    // Compute the elements of B

    B(k,k) = A(k,k)*c*c - 2.0*A(k,l)*c*s + A(l,l)*s*s;
    B(l,l) = A(l,l)*c*c + 2.0*A(k,l)*c*s + A(k,k)*s*s;

    // hard-coding non-diagonal elements
    B(k,l) = 0.0;
    B(l,k) = 0.0;

    for (int i = 0; i < N; i++)
    {
        if (i != k && i != l)
        {
            B(i,i) = A(i,i);
            B(i,k) = A(i,k)*c - A(i,l)*s;
            B(i,l) = A(i,l)*c + A(i,k)*s;
            B(k,i) = B(i,k);
            B(l,i) = B(i,l);
        }
    }



    return B;
}

// --------------------------------------------------------------------------------------------------------------
// Function for testing if the off-diagonals are small enough!


int testing(mat B, double tol)
{
    //If the off-diagonal elements are small enough,
    //testing returns 1, otherwise it returns 0

    vec max;
    max = max_element(B);
    int out;

    if (abs(max(0)) >= tol)
    {
        out = 0;
    }
    else
    {
        out = 1;
    }
    return out;
}



// -----------------------------------------------------------------------------------------------------------------

