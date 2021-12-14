---
output:
  pdf_document: default
  html_document: default
---
# DHBayes

## Purpose of Package

In settings like national-level corporate human resources, discrete data across many locations pertaining to information like site-specific employee turnover might follow a hierarchical structure.  Data like number of employees terminated in a week, or employee tenure (measured in days), might be distributed binomially or negative-binomially in a given location, respectively.  On the other hand, aggregated across many locations in a country, one can treat the parameters for each site as random variables with a common prior distribution.

This package gives tools for implementing Bayesian-oriented calculations for such discrete data, namely, maximum likelihood estimation (MLE), maximum a posteriori estimation (MAPE), empirical Bayesian (EMPB) fitting of hyperparameters, and posterior distribution sampling.  As of version 1.0.0, data which are group-level binomial, negative-binomial, and Poisson distributed are supported.  Conjugate prior forms are used for all data types.  In addition, for negative-binomially distributed data, an EMPB algorithm for a Hessian-matched Gaussian prior is provided, to across-group variability in the size parameter 'r', instead of fixed in the conjugate beta-negative-binomial scenario.

All examples in the functions' documentation generate pseudo data for testing.  As these algorithms make group-level distributional assumptions, they may not be appropriate in settings where these assumptions do not hold.  While default optimizing algorithm parameters are chosen to hopefully encourage convergence, initial starting choices can strongly influence algorithm convergence in some settings.  It is recommended to manually set starting values, where permitted, in the case where convergence fails on the default values.

## How to Install

Installation can be handled through the 'devtools' package, where installation can be done by linking to the author's github.  Copy/paste and run this code in your R console:

```{r}
install.packages('devtools')
devtools::install_github("alexandercoulter/DHBayes")
library(DHBayes)
```
\newpage
## Some Starter Example Code

Examples for each R function can be found in that respective function's documentation.  Reproduced below are the examples for the functions supporting analysis of group-level Poisson data.

### Poisson MLE

The ```mle_poisson``` function calculates the MLE of the 'lambda' parameter (L) for a collection of i.i.d. Poisson distributed data.  Copy/paste and run this code in your R console:
```{r}
# Generate example data:
set.seed(31)
L = 5

# Number of experiments, i.e. rows in df:
numexps = 10

# Filling df with pseudo data; note the requisite column 'x':
df = data.frame('x' = rpois(numexps, L))

# Generating maximum likelihood estimate (MLE) solution for L:
L_fit = mle_poisson(df = df)

# Compare fitted values to known values:
cbind(L, L_fit)
```
\newpage
### Poisson MAPE

The ```mape_poisson``` function calculates the MAPE of the 'lambda' parameter for a collection of i.i.d. Poisson distributed data, with L ~ gamma(a, b), given user-specified prior parameters 'a' and 'b'.  Copy/paste and run this code in your R console:
```{r}
# Generate example data:
set.seed(31)
L = 5

# Number of experiments, i.e. rows in df:
numexps = 10

# Filling df with pseudo data; note the requisite column 'x':
df = data.frame('x' = rpois(numexps, L))

# Generating maximum a posteriori estimate (MAPE) solution for L:
L_fit = mape_poisson(df = df, a_prior = 1, b_prior = 1)

# Compare fitted values to known values:
cbind(L, L_fit)
```
\newpage
### Poisson EMPB Prior Parameter Fitting

The ```empb_gamma_poisson_c``` function fits the prior parameters 'a' and 'b' to a collection of group-level Poisson distributed data with empirical Bayes, with the assumption that the Poisson parameters {Lg} ~ gamma(a, b).  The main numeric algorithm in this function is implemented in C++, whereas a similar function ```empb_gamma_poisson``` is entirely in R.  Copy/paste and run this code in your R console:
```{r}
# Generate example data:
set.seed(31)
a = 23
b = 9

# Number of groups:
NG = 10

# Creating group IDs:
g = replicate(NG, paste(sample(LETTERS, 10), sep="", collapse=""))

# Generating 'true' L parameters:
L = rgamma(length(g), a, b)

# Number of experiments, i.e. rows in df:
numexps = 100

# Filling df with pseudo data; note the requisite columns 'x' and 'g':
df = data.frame('x' = numeric(0), 'g' = character(0))
for(k in 1:numexps){
  gk = sample(g, 1)
  xk = rpois(1, L[g == gk])
  df = rbind(df, data.frame('x' = xk, 'g' = gk))
}

# Generating empirical Bayes (EMPB) solutions for a and b:
ab_fit = empb_gamma_poisson_c(df = df)

# Compare fitted values to known values:
cbind(c(a, b), c(ab_fit$a, ab_fit$b))
```

\newpage
### Poisson Posterior Distribution Sampling

The ```posterior_gamma_poisson``` function can generate group-level samples of the Poisson parameters {Lg} from the posterior distribution, assuming within-group data is distributed ~ poisson(Lg), and {Lg} ~ gamma(a, b).  By setting the function argument ```dist = 'pred'```, the function can instead sample from the posterior-predictive distribution.  Copy/paste and run this code in your R console:
```{r}
# Generate example data:
set.seed(31)
a = 23
b = 9

# Number of groups:
NG = 10

# Creating group IDs:
g = replicate(NG, paste(sample(LETTERS, 10), sep="", collapse=""))

# Generating 'true' L parameters:
L = rgamma(length(g), a, b)

# Number of experiments, i.e. rows in df:
numexps = 100

# Filling df with pseudo data; note the requisite columns 'x' and 'g':
df = data.frame('x' = numeric(0), 'g' = character(0))
for(k in 1:numexps){
  gk = sample(g, 1)
  xk = rpois(1, L[g == gk])
  df = rbind(df, data.frame('x' = xk, 'g' = gk))
}

# Generating 1000 posterior distribution samples for each group:
posterior_values = posterior_gamma_poisson(df = df)
dim(posterior_values)

# Create histogram of posterior distribution samples for first group (by alphabetic order):
hist(posterior_values[, 1], main = colnames(posterior_values)[1])
```
