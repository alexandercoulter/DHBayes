#' empb_beta_binomial
#'
#' @param df
#' @param eta
#' @param tol
#'
#' @return
#' @export
#'
#' @examples
empb_beta_binomial = function(df, eta = 0.1, tol = 1e-5){

  #############################################################################
  # Object 'df' should be 'data.frame', with columns 'n', 'x', and 'g'.  To that end:
  if(class(df) != 'data.frame') stop('Object \'df\' should be of type \'data.frame\'.')
  if(!all(c('n', 'x', 'g') %in% names(df))) stop('Object \'df\' must contain data.frame columns named \'n\', \'x\', and \'g\'.')
  if((length(df$n) != length(df$x)) | (length(df$n) != length(df$g))) stop('Object \'df\' data.frame columns are not of same length.')

  #############################################################################
  # Calculate parameters, starting points for optimizer:
  ug = unique(df$g)
  Lg = length(ug)
  Ng = Sg = Tg = rep(0, length(ug))
  for(j in 1:Lg){

    Sg[j] = sum(df$x[df$g == ug[j]])
    Ng[j] = sum(df$n[df$g == ug[j]])
    Tg[j] = Ng[j] - Sg[j]

  }

  #############################################################################
  # Exit if any groups have zero sample size:
  if(any(Ng == 0)) stop('Some groups have zero sample size.')

  #############################################################################
  # Initialize values:
  STg = cbind(Sg, Tg)
  ab = c(1, 1)

  #############################################################################
  # If data (df$x) are all zeros, then return(0, mean(Ng)):
  if(sum(Sg) == 0) return(c(0, mean(Ng)))
  # If data (df$x) are all ones, then return(mean(Ng), 0):
  if(sum(Tg) == 0) return(c(mean(Ng), 0))

  # Initialize parameter vectors for updating in loop:
  Score = rep(NA, 2)
  Hessian = matrix(NA, 2, 2)
  Step = c(tol, tol)

  #############################################################################
  # Run hard-capped Newton's method algorithm:
  while(sum(abs(Step)) > tol){

    #############################################################################
    # Calculate Score vector:
    Score = Lg * (digamma(sum(ab)) - digamma(ab)) + colSums(digamma(ab + STg) - digamma(sum(ab) + Ng))

    #############################################################################
    # Calculate Hessian vector:
    Hessian[ , ] = Lg * trigamma(sum(ab)) - sum(trigamma(sum(ab) + Ng))
    Hessian[c(1, 4)] = Hessian[c(1, 4)] - Lg * trigamma(ab) + colSums(trigamma(ab + STg))

    #############################################################################
    # Take damped step:
    Step = solve(Hessian, Score)
    Step = pmin(Step, 0.9 * ab)
    ab = ab - eta * Step

  }

  #############################################################################
  # Return empirical Bayes a, b:
  return(list('a' = a, 'b' = b))

}
