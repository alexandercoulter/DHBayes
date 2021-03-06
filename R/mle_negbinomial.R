#' mle_negbinomial
#'
#' @param df data.frame object, containing at least column named 'x' containing non-negative integer values.
#' @param eta positive numeric dampening parameter for Newton's method, gradient descent algorithm.
#' @param lambda non-negative numeric regularization parameter.
#' @param tol non-negative numeric tolerance parameter for exiting optimization algorithm.
#' @param maxIter positive integer setting maximum number of iterations for optimization algorithm.
#' @param method string controlling optimization method; default 'newton'.
#'
#' @return list containing MLE of negative binomial distribution parameters r and p, assuming df$x ~ nbinom(r, p).
#' @export
#'
#' @examples
#' # Generate example data:
#' set.seed(31)
#' r = 4
#' p = 0.3
#'
#' # Number of experiments, i.e. rows in df:
#' numexps = 10
#'
#' # Filling df with pseudo data; note the requisite column 'x':
#' df = data.frame('x' = rnbinom(numexps, r, p))
#'
#' # Generating maximum likelihood estimate (MLE) solution for r and p:
#' rp_fit = mle_negbinomial(df = df)
#'
#' # Compare fitted values to known values:
#' cbind(c(r, p), c(rp_fit$r, rp_fit$p))
mle_negbinomial = function(df, eta = 0.001, lambda = 0.01, tol = 0.0001, maxIter = 10000, method = c('newton', 'gdescent')){

  #############################################################################
  # Object 'df' should be 'data.frame', with column 'x'.  To that end:
  if(class(df) != 'data.frame') stop('Object \'df\' should be of type \'data.frame\'.')
  if(!('x' %in% names(df))) stop('Object \'df\' must contain data.frame column named \'x\'.')

  #############################################################################
  # Calculate reasonable starting values from MME:
  M = length(df$x)
  Sx = sum(df$x)

  if((M > 2) & (Sx > 0)){
    r0 = mean(df$x)^2 / max(var(df$x) - mean(df$x), mean(df$x))
    p0 = mean(df$x) / var(df$x)
    if(p0 >= 1) p0 = 0.9
    ell = c(log(r0), log(p0 / (1 - p0)))
  } else {
    ell = c(1, 0)
  }

  #############################################################################
  # Initialize values for Newton's method:
  Step = c(tol, tol)
  Score = rep(NA, 2)
  Hessian = matrix(NA, 2, 2)
  iternum = 0

  #############################################################################
  # Get method:
  method = match.arg(method)

  #############################################################################
  # Implement loop to calculate MLE:
  while((sum(abs(Step)) > tol) & iternum < maxIter){

    r = exp(ell[1])
    o = exp(ell[2])
    p = o / (1 + o)

    ###########################################################################
    # Cannot do Newton's method if Sx == 0 (second derivative does not exist); also do gradient descent if user requests it:
    if((Sx == 0) | (M < 4) | (method == 'gdescent')){

      #########################################################################
      # Calculate Score vector for step:
      Score = M * c(log(p), r / (o * (1 + o))) - lambda

      #########################################################################
      # Calculate step for 'ro', capping at some proportion of existing 'ro' value to not overshoot into negatives:
      Step = pmax(-0.5 * c(r, o), Score)

      #########################################################################
      # Take damped step:
      ell = log(c(r, o) + eta * Step)

    } else {

      #########################################################################
      # Calculate Score vector for Newton's step:
      Score = c(r * M * (log(p) - digamma(r)) + r * sum(digamma(df$x + r)),
                r * M * (1 - p) - p * Sx) - lambda * ell

      #########################################################################
      # Calculate Hessian Matrix for Newton's step:
      Hessian[1] = r * M * (log(p) - r * trigamma(r) - digamma(r)) + r * sum(digamma(df$'x' + r) + r * trigamma(df$'x' + r))
      Hessian[c(2, 3)] = r * (1 - p) * M
      Hessian[4] = -p * (1 - p) * (Sx + r * M)
      diag(Hessian) = diag(Hessian) - lambda

      #########################################################################
      # Take damped step:
      Step = solve(Hessian, Score)
      ell = ell - eta * Step

    }

    ###########################################################################
    # Update iteration count:
    iternum = iternum + 1

  }

  r = exp(ell[1])
  o = exp(ell[2])
  p = o / (1 + o)

  #############################################################################
  # Return MLE r, p:
  return(list('r' = r, 'p' = p))

}
