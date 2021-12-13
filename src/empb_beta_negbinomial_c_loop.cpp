// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::List empb_beta_negbinomial_c_loop(NumericVector& rab,
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
  double r = rab[0];
  double a = rab[1];
  double b = rab[2];
  NumericVector rNV(1);
  NumericVector aNV(1);
  NumericVector bNV(1);
  rNV = r;
  aNV = a;
  bNV = b;
  NumericVector R(X.size());
  NumericVector A(G);
  NumericVector B(G);

  R = r + X;
  A = a + mj * r;
  B = b + Sx;

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

    obj = sum(lbeta(A, B)) + sum(lgamma(R)) - G * (lgamma(a) + lgamma(b) - lgamma(a + b)) - Sm * lgamma(r);
    double old_obj = obj;
    double err = tol + 1;
    int iternum = 0;

    while((abs(err) > tol) && (iternum < maxIter)){

      rNV = r;
      aNV = a;
      bNV = b;

      digA   = digamma(A);
      digB   = digamma(B);
      trigA  = trigamma(A);
      trigB  = trigamma(B);
      digAB  = digamma(A + B);
      trigAB = trigamma(A + B);
      digR   = digamma(R);
      trigR  = trigamma(R);

      sum_a_digamma   = sum(digA - digAB);
      sum_b_digamma   = sum(digB - digAB);
      sum_ma_digamma  = sum(mj * (digA - digAB));
      sum_x_digamma   = sum(digR);
      digamma_gaab    = G * (digamma(aNV) - digamma(aNV + bNV))[0];
      digamma_gbab    = G * (digamma(bNV) - digamma(aNV + bNV))[0];
      sum_a_trigamma  = sum(trigA - trigAB);
      sum_b_trigamma  = sum(trigB - trigAB);
      sum_ma_trigamma = sum(mj2 * (trigA - trigAB));
      sum_x_trigamma  = sum(trigR);
      trigamma_gaab   = G * (trigamma(aNV) - trigamma(aNV + bNV))[0];
      trigamma_gbab   = G * (trigamma(bNV) - trigamma(aNV + bNV))[0];

      Score(0) = sum_ma_digamma + sum_x_digamma - Sm * digamma(rNV)[0];
      Score(1) = sum_a_digamma - digamma_gaab;
      Score(2) = sum_b_digamma - digamma_gbab;

      Hessian(0, 0) = sum_ma_trigamma + sum_x_trigamma - Sm * trigamma(rNV)[0];
      Hessian(1, 1) = sum_a_trigamma - trigamma_gaab;
      Hessian(2, 2) = sum_b_trigamma - trigamma_gbab;
      Hessian(0, 1) = sum(mj * (trigA - trigAB));
      Hessian(0, 2) = -sum(mj * trigAB);
      Hessian(1, 2) = G * trigamma(aNV + bNV)[0] - sum(trigAB);

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
      obj = sum(lbeta(A, B)) + sum(lgamma(R)) - G * (lgamma(a) + lgamma(b) - lgamma(a + b)) - Sm * lgamma(r);
      err = obj - old_obj;

      ++iternum;

    }

  } else {

    arma::colvec Score(3);
    arma::colvec Step(3);

    double sum_a_digamma;
    double sum_b_digamma;
    double sum_ma_digamma;
    double sum_x_digamma;
    double digamma_gaab;
    double digamma_gbab;

    NumericVector digA(G);
    NumericVector digB(G);
    NumericVector digAB(G);
    NumericVector digR(G);

    obj = sum(lbeta(A, B)) + sum(lgamma(R)) - G * (lgamma(a) + lgamma(b) - lgamma(a + b)) - Sm * lgamma(r);
    double old_obj = obj;
    double err = tol + 1;
    int iternum = 0;

    while((abs(err) > tol) && (iternum < maxIter)){

      rNV = r;
      aNV = a;
      bNV = b;

      digA  = digamma(A);
      digB  = digamma(B);
      digAB = digamma(A + B);
      digR  = digamma(R);

      sum_a_digamma  = sum(digA - digAB);
      sum_b_digamma  = sum(digB - digAB);
      sum_ma_digamma = sum(mj * (digA - digAB));
      sum_x_digamma  = sum(digR);
      digamma_gaab   = G * (digamma(aNV) - digamma(aNV + bNV))[0];
      digamma_gbab   = G * (digamma(bNV) - digamma(aNV + bNV))[0];

      Score(0) = sum_ma_digamma + sum_x_digamma - Sm * digamma(rNV)[0];
      Score(1) = sum_a_digamma - digamma_gaab;
      Score(2) = sum_b_digamma - digamma_gbab;

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
      obj = sum(lbeta(A, B)) + sum(lgamma(R)) - G * (lgamma(a) + lgamma(b) - lgamma(a + b)) - Sm * lgamma(r);
      err = obj - old_obj;

      ++iternum;

    }

  }

  return Rcpp::List::create(Rcpp::Named("r") = r,
                            Rcpp::Named("a") = a,
                            Rcpp::Named("b") = b,
                            Rcpp::Named("obj") = obj);

}
