#' empb_beta_negbinomial_c
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
#' # Generate example data:
#' set.seed(31)
#' r = 4
#' a = 3
#' b = 9
#'
#' # Number of groups:
#' NG = 10
#'
#' # Creating group IDs:
#' g = replicate(NG, paste(sample(LETTERS, 10), sep="", collapse=""))
#'
#' # Generating 'true' p parameters:
#' p = rbeta(length(g), a, b)
#'
#' # Number of experiments, i.e. rows in df:
#' numexps = 100
#'
#' # Filling df with pseudo data; note the requisite columns 'x' and 'g':
#' df = data.frame('x' = numeric(0), 'g' = character(0))
#' for(k in 1:numexps){
#'   gk = sample(g, 1)
#'   xk = rnbinom(1, r, p[g == gk])
#'   df = rbind(df, data.frame('x' = xk, 'g' = gk))
#' }
#'
#' # Generating empirical Bayes (EMPB) solutions for r, a, and b:
#' rab_fit = empb_beta_negbinomial_c(df = df, method = 'gdescent')
#'
#' # Compare fitted values to known values:
#' cbind(c(r, a, b), c(rab_fit$r, rab_fit$a, rab_fit$b))
empb_beta_negbinomial_c = function(df, eta = 0.1, tol = 1e-8, maxIter = 10000, starting_rab = NULL, method = c('newton', 'gdescent')){

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
  # Exit if number of groups G is 3 or fewer:
  if(G < 4) stop('Algorithm cannot provide solution for three or fewer groups.')

  #############################################################################
  # Get method:
  method = match.arg(method)

  #############################################################################
  # Give initial values for iterator:
  if(is.null(starting_rab)){
    rab = c(3, 3, 3)
  } else {
    rab = starting_rab
  }

  loop_out = empb_beta_negbinomial_c_loop(rab = rab, X = x, G = G, Sx = Sx, mj = mj, eta = eta, tol = tol, maxIter = maxIter, method = as.numeric(method != 'newton'))

  #############################################################################
  # Return empirical Bayes solutions r, a, and b:
  return(list('r' = loop_out$r, 'a' = loop_out$a, 'b' = loop_out$b))

}
