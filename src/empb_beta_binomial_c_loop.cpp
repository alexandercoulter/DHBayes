// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
#include "Rmath.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
arma::colvec empb_beta_binomial_c_loop(arma::colvec& ab,
                                       const double G,
                                       arma::rowvec& Sg,
                                       arma::rowvec& Tg,
                                       arma::rowvec& Ng,
                                       const double eta,
                                       const double tol,
                                       const int maxIter,
                                       int method) {

  if(method == 0){

    arma::mat Hessian(2, 2);

    arma::colvec Score(2);
    arma::colvec Step = ab;

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

      diab = R::digamma(sum(ab));
      triab = R::trigamma(sum(ab));
      digS = sum(digamma(as<NumericVector>(wrap(Sg + ab(0)))));
      digT = sum(digamma(as<NumericVector>(wrap(Tg + ab(1)))));
      digN = sum(digamma(as<NumericVector>(wrap(Ng + sum(ab)))));
      trigS = sum(trigamma(as<NumericVector>(wrap(Sg + ab(0)))));
      trigT = sum(trigamma(as<NumericVector>(wrap(Tg + ab(1)))));
      trigN = sum(trigamma(as<NumericVector>(wrap(Ng + sum(ab)))));

      Score(0) = G * (diab - R::digamma(ab(0))) + digS - digN;
      Score(1) = G * (diab - R::digamma(ab(1))) + digT - digN;

      Hessian(0, 0) = G * (triab - R::trigamma(ab(0))) + trigS - trigN;
      Hessian(1, 1) = G * (triab - R::trigamma(ab(1))) + trigT - trigN;
      Hessian(0, 1) = G * triab - trigN;
      Hessian(1, 0) = Hessian(0, 1);

      Step = -solve(Hessian, Score);
      ab += eta * Step;

      ++iternum;
    }

  } else {

    arma::colvec Score;
    arma::colvec Step = ab;

    double digS;
    double digT;
    double digN;
    double diab;

    int iternum = 0;

    while((sum(abs(Step)) > tol) && (iternum < maxIter)){

      diab = R::digamma(sum(ab));
      digS = sum(digamma(as<NumericVector>(wrap(Sg + ab(0)))));
      digT = sum(digamma(as<NumericVector>(wrap(Tg + ab(1)))));
      digN = sum(digamma(as<NumericVector>(wrap(Ng + sum(ab)))));

      Score(0) = G * (diab - R::digamma(ab(0))) + digS - digN;
      Score(1) = G * (diab - R::digamma(ab(1))) + digT - digN;

      Step = Score;
      ab += eta * Step;

      ++iternum;
    }

  }

  return(ab);
}
