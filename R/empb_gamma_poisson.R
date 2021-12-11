#' empb_gamma_poisson
#'
#' @param df df data.frame object, containing at least columns 'x' containing non-negative integer values, and 'g' containing group labels.
#' @param eta positive numeric dampening parameter for Newton's method, gradient descent algorithm.
#' @param tol non-negative numeric tolerance parameter for exiting optimization algorithm.
#' @param maxIter positive integer setting maximum number of iterations for optimization algorithm.
#' @param starting_ab optional 2-long numeric vector, giving initial algorithm starting point for fitting empirical Bayes estimates for 'a' and 'b'; default NULL.
#' @param method string controlling optimization method; default 'newton'.
#'
#' @return list object containing empirical Bayes (EMPB) estimates of a, b hyperparameters, assuming df$x ~ npoisson(L_g), and (L_g) ~ gamma(a, b), where 'L_g' denotes a group-level parameter.
#' @export
#'
#' @examples
empb_gamma_poisson = function(df, eta = 1, tol = 1e-10, maxIter = 10000, starting_ab = NULL, method = c('newton', 'gdescent')){

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

  #############################################################################
  # Exit if number of groups G is 2 or fewer:
  if(G < 3) stop('Algorithm cannot provide solution for three or fewer groups.')

  #############################################################################
  # Get method:
  method = match.arg(method)

  #############################################################################
  # Initialize empty objects for WHILE loop:
  Score   = rep(NA, 2)
  Hessian = matrix(NA, 2, 2)
  iternum = 0

  #############################################################################
  # Give initial values for iterator:
  if(is.null(starting_ab)){
    ab = c(mean(Sx / mj)^2 / var(Sx / mj),
           mean(Sx / mj) / var(Sx / mj))
  } else {
    if((length(starting_ab) != 2) | any(starting_ab <= 0)) stop('Vector of starting estimate \'starting_ab\' must have two positive numbers or otherwise be NULL.')
    ab = starting_ab
  }
  a = ab[1]
  b = ab[2]

  #############################################################################
  # Calculate a', b':
  a. = a + Sx
  b. = b + mj

  #############################################################################
  # Calculate initial objective function value, error:
  obj = G * (a * log(b) - lgamma(a)) + sum(lgamma(a.) - a. * log(b.))
  err = tol + 1

  #############################################################################
  # Implement while loop that fits Tau/mu by coordinate gradient/analytic descent:
  while((abs(err) > tol) & iternum < maxIter){

    a = ab[1]
    b = ab[2]

    ###########################################################################
    # Calculate a', b':
    a. = a + Sx
    b. = b + mj

    ###########################################################################
    # If/else on method input: if 'newton', perform Newton's method; if 'gdescent', perform gradient descent:
    if(method == 'newton'){

      #########################################################################
      # Calculate Score vector:
      Score = c(G * (log(b) - digamma(a)) + sum(digamma(a.) - log(b.)),
                G * a / b - sum(a. / b.))

      #########################################################################
      # Calculate Hessian matrix:
      H11 = sum(trigamma(a.)) - G * trigamma(a)
      H12 = G / b - sum(1 / b.)
      H22 = sum(a. / (b. * b.)) - G * a / (b * b)
      Hessian[1:4] = c(H11, H12, H12, H22)

      #########################################################################
      # Take damped step:
      Step = solve(Hessian, Score)
      ab = ab - eta * Step

    } else if(method == 'gdescent'){

      #########################################################################
      # Calculate Score vector:
      Score = c(G * (log(b) - digamma(a)) + sum(digamma(a.) - log(b.)),
                G * a / b - sum(a. / b.))

      #########################################################################
      # Take step:
      Step = Score
      ab = ab + eta * Step

    }

    ###########################################################################
    # Update iteration count to exit loop at maxIter:
    iternum = iternum + 1

    ###########################################################################
    # Calculate current objective function value:
    a = ab[1]
    b = ab[2]
    a. = a + Sx
    b. = b + mj
    obj_old = obj
    obj = obj = G * (a * log(b) - lgamma(a)) + sum(lgamma(a.) - a. * log(b.))

    ###########################################################################
    # Calculate change in objective function from prior step to current:
    err = obj - obj_old

  }

  #############################################################################
  # Return empirical Bayes solutions, a, and b:
  return(list('a' = a, 'b' = b, 'obj' = obj))

}


