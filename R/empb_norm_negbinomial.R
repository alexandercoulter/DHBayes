#' Title
#'
#' @param df
#' @param lambda
#' @param tol
#' @param maxIter
#' @param MLEeta
#' @param EMPBeta
#'
#' @return
#' @export
#'
#' @examples
empb_norm_negbinomial = function(df, lambda = 0.01, MLEeta = 0.1, EMPBeta = 0.01, tol = 1e-3, maxIter = 10000){
  # Object 'df' should be 'data.frame' or 'list' type, with elements 'x' and 'g'.  To that end:

  if(typeof(df) != 'list') stop('Object \'df\' should be of type \'data.frame\' or \'list\'.')
  if(!all(c('x', 'g') %in% names(df))) stop('Object \'df\' must contain list elements (or data.frame columns) named \'x\' and \'g\', respectively.')
  if(length(df$x) != length(df$g)) stop('List elements are not of same length.')

  # Dampening parameters MLEeta, EMPBeta must be positive:
  if(MLEeta <= 0) stop('Object \'MLEeta\' must be positive.')
  if(EMPBeta <= 0) stop('Object \'EMPBeta\' must be positive.')
  # Tolerance must be non-negative:
  if(tol < 0) stop('Object \'tol\' must be non-negative.')
  # Regularization parameter lambda must be non-negative:
  if(lambda < 0) stop('Object \'lambda\' must be non-negative.')
  # Maximum number of iterations must be positive integer:
  if((maxIter < 1) | ((maxIter %% 1) != 0)) stop('Object \'maxIter\' must be positive integer.')

  # Pull out unique group IDs:
  unique.g = sort(unique(df$'g'))
  G = length(unique.g)

  # Calculate mu_j's, Tau_j's for each group:
  muj = matrix(NA, nrow = G, ncol = 2)
  Tauj = array(NA, dim = c(G, 2, 2))
  mj = rep(NA, G)

  for(j in 1:G){
    d = df[df$'g' == unique.g[j], ]
    mj[j] = nrow(d)
    muj[j, ] = mle_negbinomial(df = d, eta = MLEeta, lambda = lambda, tol = tol, maxIter = maxIter)
    r = muj[j, 1]
    p = muj[j, 2]
    M = mj[j]
    Tauj[j, 1, 1] = r * M * (r * trigamma(r) + digamma(r) - log(p)) - r * sum(digamma(d$'x' + r) + r * trigamma(d$'x' + r))
    Tauj[j, 1, 2] = Tauj[j, 2, 1] = r * (p - 1) * M
    Tauj[j, 2, 2] = -1 * Tauj[j, 2, 1]
  }

  # Set muj to appropriate log, log-odds scale, instead of r/p from mle_negbinomial:
  muj[, 2] = muj[, 2] / (1 - muj[, 2])
  muj = log(muj)

  # Calculate initial values for mu, Tau:
  mu  = colSums(muj * mj) / sum(mj)
  Tau = solve(cov(muj))

  # Initialize empty objects for WHILE loop:
  Grad  = matrix(tol, 2, 2)
  Rhoji = array(NA, dim = c(G, 2, 2))
  Aj    = matrix(NA, nrow = 2, ncol = G)
  iternum = 0

  # Calculate initial objective function value:
  B = 0
  for(j in 1:G){
    Rhoji[j, , ] = solve(Tau + mj[j] * Tauj[j, , ])
    C = Tau %*% mu + mj[j] * Tauj[j, , ] %*% muj[j, ]
    B = B + t(C) %*% Rhoji[j, , ] %*% C
  }
  obj = 0.5 * (G * log(det(Tau)) - G * t(mu) %*% Tau %*% mu - sum(log(Rhoji[, 1, 1] * Rhoji[, 2, 2] - Rhoji[, 1, 2] * Rhoji[, 2, 1])) + B)
  err = tol + 1

  # Implement while loop that fits Tau/mu by coordinate gradient/analytic descent:
  while((abs(err) > tol) & iternum < maxIter){

    # Prepare calculations for new Tau gradient calculation:
    for(j in 1:G){
      Rhoji[j, , ] = solve(Tau + mj[j] * Tauj[j, , ])
      Aj[, j] = Rhoji[j, , ] %*% (Tau %*% mu + mj[j] * Tauj[j, , ] %*% muj[j, ])
    }

    # Calculate Tau gradient:
    Grad = 0.5 * (G * solve(Tau) - (colSums(Rhoji) + tcrossprod(Aj - mu)))

    # Take Tau step:
    Tau = Tau + EMPBeta * Grad

    # Calculate new mu:
    B = numeric(2)
    for(j in 1:G){
      Rhoji[j, , ] = solve(Tau + mj[j] * Tauj[j, , ])
      B = B + Rhoji[j, , ] %*% (mj[j] * Tauj[j, , ] %*% muj[j, ])
    }
    mu = c(solve(diag(2) - colSums(Rhoji) %*% Tau / G) %*% B) / G

    # Update iteration count to exit loop at maxIter:
    iternum = iternum + 1

    # Calculate current objective function value:
    B = 0
    for(j in 1:G){
      C = Tau %*% mu + mj[j] * Tauj[j, , ] %*% muj[j, ]
      B = B + t(C) %*% Rhoji[j, , ] %*% C
    }
    obj_old = obj
    obj = 0.5 * (G * log(det(Tau)) - G * t(mu) %*% Tau %*% mu - sum(log(Rhoji[, 1, 1] * Rhoji[, 2, 2] - Rhoji[, 1, 2] * Rhoji[, 2, 1])) + B)

    # Calculate change in objective function from prior step to current:
    err = obj - obj_old

  }

  return(list('mu' = mu, 'Sigma'  = solve(Tau)))
}
