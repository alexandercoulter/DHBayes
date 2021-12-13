// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::List empb_beta_negbinomial_c_loop(NumericVector rab,
                                        const NumericVector& X,
                                        const double& G,
                                        const NumericVector& Sx,
                                        const NumericVector& mj,
                                        const double& eta,
                                        const double& tol,
                                        const int maxIter,
                                        int method) {

  NumericVector mj2 = mj * mj;
  double Sm = sum(mj);
  double obj = 0;
  NumericVector r = rab[0];
  NumericVector a = rab[0];
  NumericVector b = rab[1];
  NumericVector R = r + X;
  NumericVector A = a + mj * r;
  NumericVector B = b + Sx;

  if(method == 0){

    arma::mat Hessian(3, 3);
    arma::colvec Score(3);
    arma::colvec Step(3);

    double sum_a_digamma;
    double sum_b_digamma;
    double sum_ma_digamma;
    double sum_x_digamma;
    double digamma_gaab;
    double digamma_gbab;
    double sum_a_trigamma;
    double sum_b_trigamma;
    double sum_ma_trigamma;
    double sum_x_trigamma;
    double trigamma_gaab;
    double trigamma_gbab;

    NumericVector digA(G);
    NumericVector digB(G);
    NumericVector trigA(G);
    NumericVector trigB(G);
    NumericVector digAB(G);
    NumericVector trigAB(G);
    NumericVector digR(G);
    NumericVector trigR(G);

    obj = sum(lgamma(A) + lgamma(B) - lgamma(A + B)) + sum(lgamma(R)) - G * (lgamma(a[0]) + lgamma(b[0]) - lgamma(a[0] + b[0])) - Sm * lgamma(r[0]);
    double old_obj = obj;
    double err = tol + 1;
    int iternum = 0;

    while((abs(err) > tol) && (iternum < maxIter)){

      digA = digamma(A);
      digB = digamma(B);
      trigA = trigamma(A);
      trigB = trigamma(B);
      digAB = digamma(A + B);
      trigAB = trigamma(A + B);
      digR = digamma(R);
      trigR = trigamma(R);

      sum_a_digamma = sum(digA - digAB);
      sum_b_digamma = sum(digB - digAB);
      sum_ma_digamma = sum(mj * (digA - digAB));
      sum_x_digamma = sum(digR);
      digamma_gaab = G * (digamma(a) - digamma(a + b))[0];
      digamma_gbab = G * (digamma(b) - digamma(a + b))[0];
      sum_a_trigamma = sum(trigA - trigAB);
      sum_b_trigamma = sum(trigB - trigAB);
      sum_ma_trigamma = sum(mj2 * (trigA - trigAB));
      sum_x_trigamma = sum(trigR);
      trigamma_gaab = G * (trigamma(a) - trigamma(a + b))[0];
      trigamma_gbab = G * (trigamma(b) - trigamma(a + b))[0];

      Score(0) = NULL;
      Score(1) = NULL;
      Score(2) = NULL;

      Hessian(0, 0) = NULL;
      Hessian(1, 1) = NULL;
      Hessian(2, 2) = NULL;
      Hessian(0, 1) = NULL;
      Hessian(0, 2) = NULL;
      Hessian(1, 2) = NULL;

      Hessian(1, 0) = Hessian(0, 1);
      Hessian(2, 0) = Hessian(0, 2);
      Hessian(2, 1) = Hessian(1, 2);

      Step = -solve(Hessian, Score);
      rab[0] += eta * Step(0);
      rab[1] += eta * Step(1);
      rab[2] += eta * Step(2);


      r = rab[0];
      a = rab[1];
      b = rab[2];

      R = r + X;
      A = a + mj * r;
      B = b + Sx;
      old_obj = obj;
      obj = sum(lgamma(A) + lgamma(B) - lgamma(A + B)) + sum(lgamma(R)) - G * (lgamma(a) + lgamma(b) - lgamma(a + b)) - Sm * lgamma(r);
      err = obj - old_obj;

      ++iternum;

    }

  } else {

    arma::mat Hessian(3, 3);
    arma::colvec Score(3);
    arma::colvec Step(3);

    obj = NULL;
    double old_obj = obj;
    double err = tol + 1;
    int iternum = 0;

    while((abs(err) > tol) && (iternum < maxIter)){

      Score(0) = NULL;
      Score(1) = NULL;
      Score(2) = NULL;

      Step = Score;
      rab[0] += eta * Step(0);
      rab[1] += eta * Step(1);
      rab[2] += eta * Step(2);


      r = rab[0];
      a = rab[1];
      b = rab[2];

      R = r + X;
      A = a + mj * r;
      B = b + Sx;

      old_obj = obj;
      obj = sum(lgamma(A) + lgamma(B) - lgamma(A + B)) + sum(lgamma(R)) - G * (lgamma(a) + lgamma(b) - lgamma(a + b)) - Sm * lgamma(r);
      err = obj - old_obj;

      ++iternum;

    }

  }

  return Rcpp::List::create(Rcpp::Named("r") = r,
                            Rcpp::Named("a") = a,
                            Rcpp::Named("b") = b,
                            Rcpp::Named("obj") = obj);

}
