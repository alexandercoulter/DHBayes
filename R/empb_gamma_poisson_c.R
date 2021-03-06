#' empb_gamma_poisson_c
#'
#' @param df data.frame object, containing at least columns 'x' containing non-negative integer values, and 'g' containing group labels.
#' @param eta positive numeric dampening parameter for Newton's method, gradient descent algorithm.
#' @param tol non-negative numeric tolerance parameter for exiting optimization algorithm.
#' @param maxIter positive integer setting maximum number of iterations for optimization algorithm.
#' @param starting_ab optional 2-long numeric vector, giving initial algorithm starting point for fitting empirical Bayes estimates for a and b; default NULL.
#' @param method string controlling optimization method; default 'newton'.
#'
#' @return list object containing empirical Bayes (EMPB) estimates of a, b hyperparameters, assuming df$x ~ poisson(L_g), and L_g ~ gamma(a, b), where 'L_g' denotes a group-level parameter.
#' @export
#'
#' @examples
#' # Generate example data:
#' set.seed(31)
#' a = 23
#' b = 9
#'
#' # Number of groups:
#' NG = 10
#'
#' # Creating group IDs:
#' g = replicate(NG, paste(sample(LETTERS, 10), sep="", collapse=""))
#'
#' # Generating 'true' L parameters:
#' L = rgamma(length(g), a, b)
#'
#' # Number of experiments, i.e. rows in df:
#' numexps = 100
#'
#' # Filling df with pseudo data; note the requisite columns 'x' and 'g':
#' df = data.frame('x' = numeric(0), 'g' = character(0))
#' for(k in 1:numexps){
#'   gk = sample(g, 1)
#'   xk = rpois(1, L[g == gk])
#'   df = rbind(df, data.frame('x' = xk, 'g' = gk))
#' }
#'
#' # Generating empirical Bayes (EMPB) solutions for a and b:
#' ab_fit = empb_gamma_poisson_c(df = df)
#'
#' # Compare fitted values to known values:
#' cbind(c(a, b), c(ab_fit$a, ab_fit$b))
empb_gamma_poisson_c = function(df, eta = 1, tol = 1e-8, maxIter = 10000, starting_ab = NULL, method = c('newton', 'gdescent')){

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
  x = df$x
  g = df$g

  #############################################################################
  # Set useful values for calculations, like within-group sums and sample sizes:
  for(j in 1:G){

    d = x[g == unique.g[j]]
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
  # Give initial values for iterator:
  if(is.null(starting_ab)){
    msx = mean(Sx)
    vsx = var(Sx)
    ab = pmax(c(msx * msx / vsx, msx / vsx), 2)
  } else {
    if((length(starting_ab) != 2) | any(starting_ab <= 0)) stop('Vector of starting estimate \'starting_ab\' must have two positive numbers or otherwise be NULL.')
    ab = starting_ab
  }

  loop_out = empb_gamma_poisson_c_loop(ab = ab, G = G, Sx = Sx, mj = mj, eta = eta, tol = tol, maxIter = maxIter, method = as.numeric(method != 'newton'))

  #############################################################################
  # Return empirical Bayes solutions, a, and b:
  return(list('a' = loop_out$a,
              'b' = loop_out$b,
              'obj' = loop_out$obj))

}


