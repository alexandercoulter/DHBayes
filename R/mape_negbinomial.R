#' Title
#'
#' @param df
#' @param r_prior
#' @param p_prior
#' @param M_prior
#' @param eta
#' @param tol
#' @param maxIter
#'
#' @return
#' @export
#'
#' @examples
mape_negbinomial = function(df, r_prior, p_prior, M_prior, eta = 1, tol = 0.0001, maxIter = 200){

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
    # Calculate Score vector for Newton's step:
    Score = c(r * M * (log(p) - digamma(r)) + r * sum(digamma(df$x + r)),
              r * M * (1 - p) - p * Sx)
    Score = Score - Sigma_i %*% (ell - mu)

    ###########################################################################
    # Calculate Hessian Matrix for Newton's step:
    Hessian[1] = r * M * (log(p) - digamma(r) - r * trigamma(r)) + r * sum(digamma(df$x + r) + r * trigamma(df$x + r))
    Hessian[c(2, 3)] = r * (1 - p) * M
    Hessian[4] = -1 * p * (1 - p) * (Sx + r * M)
    Hessian = Hessian - Sigma_i

    ###########################################################################
    # Calculate step for 'ro':
    Step = solve(Hessian, Score)

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
  return(c(r, p))

}
