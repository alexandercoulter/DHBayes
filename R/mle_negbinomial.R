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
mle_negbinomial = function(df, eta = 1, lambda = 0.01, tol = 0.0001, maxIter = 200){
  # Object 'df' should be 'data.frame' or 'list' type, with element 'x'.  To that end:

  if(typeof(df) != 'list') stop('Object \'df\' should be of type \'data.frame\' or \'list\'.')
  if(!('x' %in% names(df))) stop('Object \'df\' must contain list element (or data.frame column) named \'x\'.')

  # Calculate MLE from input:
  M = length(df$x)
  Sx = sum(df$x)
  ro = c(1, M / Sx)

  Step = c(tol, tol)
  Score = rep(NA, 2)
  Hessian = matrix(NA, 2, 2)
  iternum = 0

  while((sum(abs(Step)) > tol) & iternum < maxIter){
    r = ro[1]
    o = ro[2]

    # Calculate Score vector for Newton's step:
    Score = c(M * (log(o / (1 + o)) - digamma(r)) + sum(digamma(df$x + r)) - lambda,
              sum(r / o - (df$x + r) / (1 + o)))

    # Calculate Hessian Matrix for Newton's step:
    Hessian[1] = sum(trigamma(df$x + r)) - M * trigamma(r)
    Hessian[c(2, 3)] = M / (o * (1 + o))
    Hessian[4] = sum((df$x + r) / (1 + o)^2 - r / o^2)

    # Calculate step for 'ro', capping at some proportion of existing 'ro' value to not overshoot into negatives:
    Step = pmin(0.5 * ro, solve(Hessian, Score))

    # Take damped step:
    ro = ro - eta * Step

    iternum = iternum + 1
  }
  return(c(ro[1], ro[2] / (1 + ro[2])))
}
