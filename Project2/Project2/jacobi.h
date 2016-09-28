#ifndef JACOBI_H
#define JACOBI_H
#include <armadillo>

arma::vec max_element(arma::mat A);
arma::vec trigonometry(arma::mat A, int k, int l);
arma::mat make_B(arma::mat A, int k, int l, double c, double s);
int testing(arma::mat B, double tol=1e-8);




#endif // JACOBI_H
