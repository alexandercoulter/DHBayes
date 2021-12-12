// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector empb_beta_binomial_c_loop(const double& G,
                                        const NumericVector& Sg,
                                        const NumericVector& Tg,
                                        const NumericVector& Ng,
                                        const double& eta,
                                        const double& tol,
                                        const int maxIter,
                                        int method) {

  NumericVector ab(2);
  ab[0] = 1;
  ab[1] = 1;

  if(method == 0){

    arma::mat Hessian(2, 2);
    NumericVector digaORb(2);
    NumericVector trigaORb(2);

    arma::colvec Score(2);
    arma::colvec Step(2);
    Step(0) = ab[0];
    Step(1) = ab[1];

    double digS;
    double digT;
    double digN;
    double trigS;
    double trigT;
    double trigN;
    double diab;
    double triab;

    int iternum = 0;

    while((sum(abs(Step)) > tol) && (iternum < maxIter)){

      digaORb = digamma(ab);
      trigaORb = trigamma(ab);
      diab = digamma(NumericVector(1, sum(ab)))[0];
      triab = trigamma(NumericVector(1, sum(ab)))[0];
      digS = sum(digamma(Sg + ab[0]));
      digT = sum(digamma(Tg + ab[1]));
      digN = sum(digamma(Ng + sum(ab)));
      trigS = sum(trigamma(Sg + ab[0]));
      trigT = sum(trigamma(Tg + ab[1]));
      trigN = sum(trigamma(Ng + sum(ab)));

      Score(0) = G * (diab - digaORb[0]) + digS - digN;
      Score(1) = G * (diab - digaORb[1]) + digT - digN;

      Hessian(0, 0) = G * (triab - trigaORb[0]) + trigS - trigN;
      Hessian(1, 1) = G * (triab - trigaORb[1]) + trigT - trigN;
      Hessian(0, 1) = G * triab - trigN;
      Hessian(1, 0) = Hessian(0, 1);

      Step = -solve(Hessian, Score);
      ab[0] += eta * Step(0);
      ab[1] += eta * Step(1);

      ++iternum;
    }

  } else {

    NumericVector digaORb(2);
    NumericVector trigaORb(2);

    arma::colvec Score(2);
    arma::colvec Step(2);
    Step(0) = ab[0];
    Step(1) = ab[1];

    double digS;
    double digT;
    double digN;
    double diab;

    int iternum = 0;

    while((sum(abs(Step)) > tol) && (iternum < maxIter)){

      digaORb = digamma(ab);
      diab = digamma(NumericVector(1, sum(ab)))[0];
      digS = sum(digamma(Sg + ab[0]);
      digT = sum(digamma(Tg + ab[1]);
      digN = sum(digamma(Ng + sum(ab)));

      Score(0) = G * (diab - digaORb[0]) + digS - digN;
      Score(1) = G * (diab - digaORb[1]) + digT - digN;

      Step = Score;
      ab[0] += eta * Step(0);
      ab[1] += eta * Step(1);

      ++iternum;
    }

  }

  return(ab);
}
