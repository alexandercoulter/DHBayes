#' posterior_beta_binomial
#'
#' @param df data.frame object, containing at least columns named 'x' containing non-negative integer values (number of successes), 'n' containing non-negative integer values (number of trials), and 'g' containing group labels.
#' @param ab_prior
#' @param dist
#' @param Nsamp
#' @param eta positive numeric dampening parameter for Newton's method, gradient descent algorithm.
#' @param tol non-negative numeric tolerance parameter for exiting optimization algorithm.
#' @param maxIter positive integer setting maximum number of iterations for optimization algorithm.
#' @param method string controlling optimization method; default 'newton'.

#'
#' @return list object containing empirical Bayes (EMPB) estimates of a, b hyperparameters, assuming df$x ~ binom(p_g, df$n), and p_g ~ beta(a, b), where 'p_g' denotes a group-level parameter.
#' @export
#'
#' @examples
#' # Hello world
#' x = 2
posterior_beta_binomial = function(df, ab_prior = NULL, dist = c('post', 'pred'), Nsamp = 1000, ...){

  #############################################################################
  # Object 'df' should be 'data.frame', with columns 'n', 'x', and 'g'.  To that end:
  if(class(df) != 'data.frame') stop('Object \'df\' should be of type \'data.frame\'.')
  if(!all(c('n', 'x', 'g') %in% names(df))) stop('Object \'df\' must contain data.frame columns named \'n\', \'x\', and \'g\'.')
  if((length(df$n) != length(df$x)) | (length(df$n) != length(df$g))) stop('Object \'df\' data.frame columns are not of same length.')

  #############################################################################
  # If 'ab_prior' is NULL, then calculate prior a, b parameters from empirical Bayes; otherwise, extract them:
  if(is.null(ab_prior)){
    ab = empb_beta_binomial(df = df, ...)
    a = ab$a
    b = ab$b
  } else {
    a = ab_prior[1]
    b = ab_prior[2]
  }

  #############################################################################
  # Calculate posterior parameter values:
  a. = a + sum(df$x)
  b. = b + sum(df$n) - (a. - a)

  #############################################################################
  # If 'dist' is 'post' (for 'posterior distribution'), then generate samples from posterior distribution; otherwise, if 'dist' is 'pred' (for 'posterior predictive distribution'), then generate samples from that:
  dist = match.arg(dist)
  if(dist == 'post'){

    ###########################################################################
    # Generate samples:
    X = NULL

  } else {

    ###########################################################################
    # Generate samples:
    X = NULL

  }

  ###########################################################################
  # Set names of output matrix:
  names(X) = sort(unique(df$'g'))

  #############################################################################
  # Return samples:
  return(X)

}


X = matrix(1:9, 3, 3)
colnames(X) = letters[1:3]
