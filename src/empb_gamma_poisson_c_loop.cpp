// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector empb_gamma_poisson_c_loop(NumericVector ab,
                                        const double& G,
                                        const NumericVector& Sx,
                                        const NumericVector& mj
                                        const double& eta,
                                        const double& tol,
                                        const int maxIter,
                                        int method) {

  if(method == 0){

    arma::mat Hessian(2, 2);
    NumericVector digaORb(2);
    NumericVector trigaORb(2);

    arma::colvec Score(2);
    arma::colvec Step(2);
    double a = ab[0];
    double b = ab[1];
    NumericVector A = a + Sx;
    NumericVector B = b + mj;

    double obj = G * (a * log(b) - lgamma(a)) + sum(lgamma(A) - A * log(B));
    double err = tol + 1;
    int iternum = 0;

    while((abs(err) > tol) && (iternum < maxIter)){

      a = ab[0];
      b = ab[1];
      A = a + Sx;
      B = b + mj;

      Score(0) = G * (log(b) - digamma(NumericVector(1, a))[0]) + sum(digamma(A) - log(B));
      Score(1) = G * a / b - sum(A / B);

      Hessian(0, 0) =
      Hessian(1, 1) =
      Hessian(0, 1) =
      Hessian(1, 0) =

      Step = -solve(Hessian, Score);
      ab[0] += eta * Step(0);
      ab[1] += eta * Step(1);

      ++iternum;

    }

  } else {

    arma::mat Hessian(2, 2);
    NumericVector digaORb(2);
    NumericVector trigaORb(2);

    arma::colvec Score(2);
    arma::colvec Step(2);
    double a = ab[0];
    double b = ab[1];
    NumericVector A = a + Sx;
    NumericVector B = b + mj;

    double obj = G * (a * log(b) - lgamma(a)) + sum(lgamma(A) - A * log(B));
    double err = tol + 1;
    int iternum = 0;

    while((abs(err) > tol) && (iternum < maxIter)){

      a = ab[0];
      b = ab[1];
      A = a + Sx;
      B = b + mj;

      Score(0) = G * (log(b) - digamma(NumericVector(1, a))[0]) + sum(digamma(A) - log(B));
      Score(1) = G * a / b - sum(A / B);

      Step = Score;
      ab[0] += eta * Step(0);
      ab[1] += eta * Step(1);

      ++iternum;

    }

  }

  return(ab);
}
