#' empb_beta_negbinomial
#'
#' @param df data.frame object, containing at least columns 'x' containing non-negative integer values, and 'g' containing group labels.
#' @param eta positive numeric dampening parameter for Newton's method, gradient descent algorithm.
#' @param tol non-negative numeric tolerance parameter for exiting optimization algorithm.
#' @param maxIter positive integer setting maximum number of iterations for optimization algorithm.
#' @param starting_rab optional 3-long numeric vector, giving initial algorithm starting point for fitting empirical Bayes estimates for r, a, and b, respectively; default NULL.
#' @param method string controlling optimization method; default 'newton'.
#'
#' @return list object containing empirical Bayes (EMPB) estimates of r, a, and b hyperparameters, assuming df$x ~ nbinom(r, p_g), and p_g ~ beta(a, b), where 'p_g' denotes a group-level parameter.
#' @export
#'
#' @examples
empb_beta_negbinomial = function(df, eta = 1, tol = 1e-8, maxIter = 10000, starting_rab = NULL, method = c('newton', 'gdescent')){

  #############################################################################
  # Object 'df' should be 'data.frame', with columns 'x' and 'g'.  To that end:
  if(class(df) != 'data.frame') stop('Object \'df\' should be of type \'data.frame\'.')
  if(!all(c('x', 'g') %in% names(df))) stop('Object \'df\' must contain data.frame columns named \'x\' and \'g\'.')
  if(length(df$x) != length(df$g)) stop('Object \'df\' data.frame columns \'x\' and \'g\' are not of same length.')

  #############################################################################
  # Dampening parameter eta must be positive:
  if(eta <= 0) stop('Object \'eta\' must be positive.')
  # Tolerance must be non-negative:
  if(tol < 0) stop('Object \'tol\' must be non-negative.')
  # Maximum number of iterations must be positive integer:
  if((maxIter < 1) | ((maxIter %% 1) != 0)) stop('Object \'maxIter\' must be positive integer.')

  #############################################################################
  # Pull out unique group IDs:
  unique.g = sort(unique(df$'g'))
  G = length(unique.g)
  Sx = mj = rep(NA, G)
  x = df$'x'

  #############################################################################
  # Set useful values for calculations, like within-group sums and sample sizes:
  for(j in 1:G){

    d = x[df$'g' == unique.g[j]]
    Sx[j] = sum(d)
    mj[j] = length(d)

  }
  Sm = sum(mj)

  #############################################################################
  # Exit if number of groups G is 3 or fewer:
  if(G < 4) stop('Algorithm cannot provide solution for three or fewer groups.')

  #############################################################################
  # Get method:
  method = match.arg(method)

  #############################################################################
  # Initialize empty objects for WHILE loop:
  Score  = rep(NA, 3)
  Hessian = matrix(NA, 3, 3)
  iternum = 0

  #############################################################################
  # Give initial values for iterator:
  if(is.null(starting_rab)){
    rab = c(1, 1, 1)
  } else {
    rab = starting_rab
  }
  r = rab[1]
  a = rab[2]
  b = rab[3]

  #############################################################################
  # Calculate a', b':
  a. = a + mj * r
  b. = b + Sx

  #############################################################################
  # Calculate initial objective function value, error:
  obj = sum(lgamma(x + r)) + sum(lbeta(a., b.)) - G * lbeta(a, b) - Sm * lgamma(r)
  err = tol + 1

  #############################################################################
  # Implement while loop that fits Tau/mu by coordinate gradient/analytic descent:
  while((abs(err) > tol) & iternum < maxIter){

    r = rab[1]
    a = rab[2]
    b = rab[3]

    ###########################################################################
    # Calculate a', b':
    a. = a + mj * r
    b. = b + Sx

    ###########################################################################
    # If/else on method input: if 'newton', perform Newton's method; if 'gdescent', perform gradient descent:
    if(method == 'newton'){

      #########################################################################
      # Calculate preparatory values:
      sum_a_digamma = sum(digamma(a.) - digamma(a. + b.))
      sum_b_digamma = sum(digamma(b.) - digamma(a. + b.))
      sum_ma_digamma = sum(mj * (digamma(a.) - digamma(a. + b.)))
      sum_x_digamma = sum(digamma(x + r))
      digamma_gaab = G * (digamma(a) - digamma(a + b))
      digamma_gbab = G * (digamma(b) - digamma(a + b))
      sum_a_trigamma = sum(trigamma(a.) - trigamma(a. + b.))
      sum_b_trigamma = sum(trigamma(b.) - trigamma(a. + b.))
      sum_ma_trigamma = sum(mj^2 * (trigamma(a.) - trigamma(a. + b.)))
      sum_x_trigamma = sum(trigamma(x + r))
      trigamma_gaab = G * (trigamma(a) - trigamma(a + b))
      trigamma_gbab = G * (trigamma(b) - trigamma(a + b))

      #########################################################################
      # Calculate Score vector:
      Score = c(sum_ma_digamma + sum_x_digamma - Sm * digamma(r),
                sum_a_digamma - digamma_gaab,
                sum_b_digamma - digamma_gbab)

      #########################################################################
      # Calculate Hessian components:
      H11 = sum_ma_trigamma + sum_x_trigamma - Sm * trigamma(r)
      H22 = sum_a_trigamma - trigamma_gaab
      H33 = sum_b_trigamma - trigamma_gbab
      H12 = sum(mj * (trigamma(a.) - trigamma(a. + b.)))
      H13 = -sum(mj * trigamma(a. + b.))
      H23 = G * trigamma(a + b) - sum(trigamma(a. + b.))

      #########################################################################
      # Combine Hessian components into Hessian matrix:
      Hessian = matrix(c(H11, H12, H13, H12, H22, H23, H13, H23, H33), 3, 3)

      #########################################################################
      # Take damped step:
      Step = solve(Hessian, Score)
      rab = rab - eta * Step

    } else if(method == 'gdescent'){

      #########################################################################
      # Calculate preparatory values:
      sum_a_digamma = sum(digamma(a.) - digamma(a. + b.))
      sum_b_digamma = sum(digamma(b.) - digamma(a. + b.))
      sum_ma_digamma = sum(mj * (digamma(a.) - digamma(a. + b.)))
      sum_x_digamma = sum(digamma(x + r))
      digamma_gaab = G * (digamma(a) - digamma(a + b))
      digamma_gbab = G * (digamma(b) - digamma(a + b))

      #########################################################################
      # Calculate Score vector:
      Score = c(sum_ma_digamma + sum_x_digamma - Sm * digamma(r),
                sum_a_digamma - digamma_gaab,
                sum_b_digamma - digamma_gbab)

      #########################################################################
      # Take step:
      Step = Score
      rab = rab + eta * Step

    }

    ###########################################################################
    # Update iteration count to exit loop at maxIter:
    iternum = iternum + 1

    ###########################################################################
    # Calculate current objective function value:
    r = rab[1]
    a = rab[2]
    b = rab[3]
    a. = a + mj * r
    b. = b + Sx
    obj_old = obj
    obj = sum(lgamma(x + r)) + sum(lbeta(a., b.)) - G * lbeta(a, b) - Sm * lgamma(r)

    ###########################################################################
    # Calculate change in objective function from prior step to current:
    err = obj - obj_old

  }

  #############################################################################
  # Return empirical Bayes solutions r, a, and b:
  return(list('r' = r, 'a' = a, 'b' = b))

}
