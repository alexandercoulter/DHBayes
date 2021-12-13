// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::List empb_gamma_poisson_c_loop(NumericVector ab,
                                     const double& G,
                                     const NumericVector& Sx,
                                     const NumericVector& mj,
                                     const double& eta,
                                     const double& tol,
                                     const int maxIter,
                                     int method) {

  double obj = 0;
  arma::colvec Score(2);
  arma::colvec Step(2);
  double a = ab[0];
  double b = ab[1];
  NumericVector aNV(1);
  aNV = a;
  NumericVector A = a + Sx;
  NumericVector B = b + mj;

  if(method == 0){

    arma::mat Hessian(2, 2);

    obj = G * (a * log(b) - lgamma(a)) + sum(lgamma(A) - A * log(B));
    double old_obj = obj;
    double err = tol + 1;
    int iternum = 0;

    while((abs(err) > tol) && (iternum < maxIter)){

      aNV = a;

      Score(0) = G * (log(b) - digamma(aNV)[0]) + sum(digamma(A) - log(B));
      Score(1) = G * a / b - sum(A / B);

      Hessian(0, 0) = sum(trigamma(A)) - G * trigamma(aNV)[0];
      Hessian(1, 1) = sum(A / (B * B)) - G * a / (b * b);
      Hessian(0, 1) = G / b - sum(1 / B);
      Hessian(1, 0) = Hessian(0, 1);

      Step = -solve(Hessian, Score);
      ab[0] += eta * Step(0);
      ab[1] += eta * Step(1);

      a = ab[0];
      b = ab[1];
      A = a + Sx;
      B = b + mj;
      old_obj = obj;
      obj = G * (a * log(b) - lgamma(a)) + sum(lgamma(A) - A * log(B));
      err = obj - old_obj;

      ++iternum;

    }

  } else {

    obj = G * (a * log(b) - lgamma(a)) + sum(lgamma(A) - A * log(B));
    double old_obj = obj;
    double err = tol + 1;
    int iternum = 0;

    while((abs(err) > tol) && (iternum < maxIter)){

      aNV = a;

      Score(0) = G * (log(b) - digamma(aNV)[0]) + sum(digamma(A) - log(B));
      Score(1) = G * a / b - sum(A / B);

      Step = Score;
      ab[0] += eta * Step(0);
      ab[1] += eta * Step(1);

      a = ab[0];
      b = ab[1];
      A = a + Sx;
      B = b + mj;
      old_obj = obj;
      obj = G * (a * log(b) - lgamma(a)) + sum(lgamma(A) - A * log(B));
      err = obj - old_obj;

      ++iternum;

    }

  }

  return Rcpp::List::create(Rcpp::Named("a") = a,
                            Rcpp::Named("b") = b,
                            Rcpp::Named("obj") = obj);

}
