---
output:
  pdf_document: default
  html_document: default
---
# DHBayes

### Updates Through 23 November 2021

The "DHBayes" package currently implements maximum likelihood estimation (MLE) for binomially distributed data X_i ~ binom(n_i, p), maximum a posteriori estimation (MAPE) for such binomially distributed data with a (conjugate) prior ~ beta(a, b) on 'p'; empirical Bayesian fitting of the prior parameter 'a' and 'b' in the previous conjugate prior model for binomially distributed data.  In this latter scenario, the data are required to be indexed by a group ID, which are permitted to be repeated among different observations, fitting an independent-group hierarchical setting.  These functions have compatibility checks, and are planned to be publicy available to the user upon download.

### To Finish

These same functionalities—MLE, MAPE for a conjugate prior model, and empirical Bayesian fitting for prior parameters in a hierarchical Bayesian setting—will be implemented for negative-binomially distributed data (with a beta conjugate prior) and Poisson distributed data (with a gamma conjugate prior).  Further, for a fully Bayesian model, hyperprior parameter estimation for conjugate priors for these prior distributions' parameters will be available, again through empirical Bayesian methods.  The principal function this package will provide, however, will be sampling algorithms for these hierarchical Bayesian models.  These samplers will support options for user-provided hyperprior parameter values, or specifying which parameters should be fit by empirical Bayes while the others are run through a sampler.

Numerical estimation algorithms, like many of the MLE, MAPE, and empirical Bayes estimation functions, will ultimately be implemented in C++.  The sampler algorithms will also be implemented in C++.
