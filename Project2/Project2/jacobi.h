#ifndef JACOBI_H
#define JACOBI_H
#include <armadillo>

arma::vec max_element(arma::mat A);
arma::vec trigonometry(arma::mat A, int k, int l);
arma::mat make_B(arma::mat A, int k, int l, double c, double s);


#endif // JACOBI_H
