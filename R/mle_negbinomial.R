#' Title
#'
#' @param df
#' @param eta
#' @param lambda
#' @param tol
#' @param maxIter
#'
#' @return
#' @export
#'
#' @examples
mle_negbinomial = function(df, eta = 0.001, lambda = 0.01, tol = 0.0001, maxIter = 10000){
  # Object 'df' should be 'data.frame' or 'list' type, with element 'x'.  To that end:

  if(typeof(df) != 'list') stop('Object \'df\' should be of type \'data.frame\' or \'list\'.')
  if(!('x' %in% names(df))) stop('Object \'df\' must contain list element (or data.frame column) named \'x\'.')

  # Calculate MLE from input:
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

  Step = c(tol, tol)
  Score = rep(NA, 2)
  Hessian = matrix(NA, 2, 2)
  iternum = 0

  while((sum(abs(Step)) > tol) & iternum < maxIter){
    r = exp(ell[1])
    o = exp(ell[2])
    p = o / (1 + o)

    # Cannot do Newton's method if Sx == 0 (second derivative does not exist):
    if((Sx == 0) | (M < 4)){

      # Calculate Score vector for step:
      Score = M * c(log(p), r / (o * (1 + o))) - lambda

      # Calculate step for 'ro', capping at some proportion of existing 'ro' value to not overshoot into negatives:
      Step = pmax(-0.5 * c(r, o), Score)

      # Take damped step:
      ell = log(c(r, o) + eta * Step)

    } else {

      # Calculate Score vector for Newton's step:
      Score = c(r * M * (log(p) - digamma(r)) + r * sum(digamma(df$x + r)),
                r * M * (1 - p) - p * Sx) - lambda * ell

      # Calculate Hessian Matrix for Newton's step:
      Hessian[1] = r * M * (log(p) - r * trigamma(r) - digamma(r)) + r * sum(digamma(df$'x' + r) + r * trigamma(df$'x' + r))
      Hessian[c(2, 3)] = r * (1 - p) * M
      Hessian[4] = -p * (1 - p) * (Sx + r * M)
      diag(Hessian) = diag(Hessian) - lambda

      # Take damped step:
      Step = solve(Hessian, Score)
      ell = ell - eta * Step

    }

    # Update iteration count:
    iternum = iternum + 1

  }

  r = exp(ell[1])
  o = exp(ell[2])
  p = o / (1 + o)

  return(c(r, p))
}
