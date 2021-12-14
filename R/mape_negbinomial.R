#' mape_negbinomial
#'
#' @param df data.frame object, containing at least column named 'x' containing non-negative integer values.
#' @param r_prior positive numeric; recommended choose as if there was observed "prior sample" Y of size 'M_prior', where Y ~ nbinom(r_prior, p_prior).
#' @param p_prior numeric strictly greater than 0 and less than 1; recommended choose as if there was observed "prior sample" Y of size 'M_prior', where Y ~ nbinom(r_prior, p_prior).
#' @param M_prior positive numeric; recommended choose as if there was observed "prior sample" Y of size 'M_prior', where Y ~ nbinom(r_prior, p_prior).
#' @param eta positive numeric dampening parameter for Newton's method, gradient descent algorithm.
#' @param tol non-negative numeric tolerance parameter for exiting optimization algorithm.
#' @param maxIter positive integer setting maximum number of iterations for optimization algorithm.
#' @param method string controlling optimization method; default 'newton'.
#'
#' @return List object containing maximum a-posteriori estimates (MAPEs) of negative binomial distribution parameters r and p, assuming df$x ~ nbinom(r, p), and (r, p) ~ MVN(mu, Sigma).  (Please see DHBayes_Derivations.pdf on GitHub for how mu, Sigma relate to 'r_prior', 'p_prior', and 'M_prior' inputs.)
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
#' # Generating maximum a posteriori estimate (MAPE) solution for p:
#' rp_fit = mape_negbinomial(df = df, r_prior = 2, p_prior = 0.5, M_prior = 1)
#'
#' # Compare fitted values to known values:
#' cbind(c(r, p), c(rp_fit$r, rp_fit$p))
mape_negbinomial = function(df, r_prior, p_prior, M_prior, eta = 0.1, tol = 0.0001, maxIter = 10000, method = c('newton', 'gdescent')){

  #############################################################################
  # Object 'df' should be 'data.frame', with column 'x'.  To that end:
  if(class(df) != 'data.frame') stop('Object \'df\' should be of type \'data.frame\'.')
  if(!('x' %in% names(df))) stop('Object \'df\' must contain data.frame column named \'x\'.')

  #############################################################################
  # Compatibility checks:
  # Prior parameter r_prior must be positive:
  if(r_prior <= 0) stop('Object \'r_prior\' must be positive.')
  # Prior parameter p_prior must be strictly between 0 and 1:
  if((p_prior <= 0) | (p_prior >= 1)) stop('Object \'p_prior\' must be strictly between zero and one.')
  # Prior parameter M_prior must be positive:
  if(M_prior <= 0) stop('Object \'M_prior\' must be positive.')
  # Dampening parameter eta must be positive:
  if(eta <= 0) stop('Object \'eta\' must be positive.')
  # Tolerance must be non-negative:
  if(tol < 0) stop('Object \'tol\' must be non-negative.')
  # Maximum number of iterations must be positive integer:
  if((maxIter < 1) | ((maxIter %% 1) != 0)) stop('Object \'maxIter\' must be positive integer.')

  #############################################################################
  # Set prior mean vector mu, starting point:
  mu = log(c(r_prior, p_prior / (1 - p_prior)))
  ell = mu

  #############################################################################
  # Calculate prior precision matrix Sigma^(-1):
  h11 = r_prior * M_prior * (sum(dnbinom(0:qnbinom(0.99999, r_prior, p_prior), r_prior, p_prior) * digamma(0:qnbinom(0.99999, r_prior, p_prior) + r_prior)) - digamma(r_prior) + log(p_prior) +
                              r * (sum(dnbinom(0:qnbinom(0.99999, r_prior, p_prior), r_prior, p_prior) * trigamma(0:qnbinom(0.99999, r_prior, p_prior) + r_prior)) - trigamma(r_prior)))
  h234 = r_prior * (1 - p_prior) * M_prior
  Sigma_i = matrix(c(-h11, -h234, -h234, h234), 2, 2)

  #############################################################################
  # Calculate some values for score, Hessian:
  M = length(df$x)
  Sx = sum(df$x)

  #############################################################################
  # Get method:
  method = match.arg(method)

  #############################################################################
  # Initialize empty objects for WHILE loop:
  Step = c(tol, tol)
  Score = rep(NA, 2)
  Hessian = matrix(NA, 2, 2)
  iternum = 0

  #############################################################################
  # Implement loop to solve MAPE:
  while((sum(abs(Step)) > tol) & iternum < maxIter){

    r = exp(ell[1])
    o = exp(ell[2])
    p = o / (1 + o)

    ###########################################################################
    # Calculate Score vector:
    Score = c(r * M * (log(p) - digamma(r)) + r * sum(digamma(df$x + r)),
              r * M * (1 - p) - p * Sx)
    Score = Score - Sigma_i %*% (ell - mu)

    ###########################################################################
    # Perform Newton's method step if method == 'newton', otherwise gradient descent:
    if(method == 'newton'){

      #########################################################################
      # Calculate Hessian Matrix for Newton's step if method = 'newton':
      Hessian[1] = r * M * (log(p) - digamma(r) - r * trigamma(r)) + r * sum(digamma(df$x + r) + r * trigamma(df$x + r))
      Hessian[c(2, 3)] = r * (1 - p) * M
      Hessian[4] = -1 * p * (1 - p) * (Sx + r * M)
      Hessian = Hessian - Sigma_i

      #########################################################################
      # Calculate step for 'ro':
      Step = solve(Hessian, Score)

    } else {

      #########################################################################
      # Calculate step for 'ro':
      Step = -Score

    }

    ###########################################################################
    # Take damped step:
    ell = ell - eta * Step

    ###########################################################################
    # Update iteration count to allow early exit:
    iternum = iternum + 1

  }

  r = exp(ell[1])
  o = exp(ell[2])
  p = o / (1 + o)

  #############################################################################
  # Return MAPE r, p:
  return(list('r' = r, 'p' = p))

}
