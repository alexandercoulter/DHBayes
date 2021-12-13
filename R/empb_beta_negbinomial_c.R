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
#' # Hello world
#' x = 2
empb_beta_negbinomial_c = function(df, eta = 1, tol = 1e-8, maxIter = 10000, starting_rab = NULL, method = c('newton', 'gdescent')){

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

  loop_out = empb_beta_negbinomial_c_loop(rab = rab, X = x, G = G, Sx = Sx, mj = mj, eta = eta, tol = tol, maxIter = maxIter, method = as.numeric(method != 'newton'))

  #############################################################################
  # Return empirical Bayes solutions r, a, and b:
  return(list('r' = loop_out$r, 'a' = loop_out$a, 'b' = loop_out$b))

}
