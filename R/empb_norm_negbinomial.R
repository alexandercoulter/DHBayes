#' empb_norm_negbinomial
#'
#' @param df
#' @param lambda
#' @param MLEeta
#' @param EMPBeta
#' @param tol
#' @param maxIter
#' @param jitter
#'
#' @return
#' @export
#'
#' @examples
empb_norm_negbinomial = function(df, lambda = 1, MLEeta = 0.001, EMPBeta = 0.001, tol = 0.1, maxIter = 10000, jitter = FALSE){

  #############################################################################
  # Object 'df' should be 'data.frame', with columns 'x' and 'g'.  To that end:
  if(class(df) != 'data.frame') stop('Object \'df\' should be of type \'data.frame\'.')
  if(!all(c('x', 'g') %in% names(df))) stop('Object \'df\' must contain data.frame columns named \'x\' and \'g\'.')
  if(length(df$x) != length(df$g)) stop('Object \'df\' data.frame columns \'x\' and \'g\' are not of same length.')

  #############################################################################
  # Dampening parameters MLEeta, EMPBeta must be positive:
  if(MLEeta <= 0) stop('Object \'MLEeta\' must be positive.')
  if(EMPBeta <= 0) stop('Object \'EMPBeta\' must be positive.')
  # Tolerance must be non-negative:
  if(tol < 0) stop('Object \'tol\' must be non-negative.')
  # Regularization parameter lambda must be non-negative:
  if(lambda < 0) stop('Object \'lambda\' must be non-negative.')
  # Maximum number of iterations must be positive integer:
  if((maxIter < 1) | ((maxIter %% 1) != 0)) stop('Object \'maxIter\' must be positive integer.')

  #############################################################################
  # Pull out unique group IDs:
  unique.g = sort(unique(df$'g'))
  G = length(unique.g)

  #############################################################################
  # Exit if number of groups G is 1; exit and provide tip for forcing solution is number of groups G is 2:
  if(G == 1) stop('Algorithm cannot provide solution for one group.')
  if(G == 2) print('WARNING: algorithm result very unreliable for two groups.')

  #############################################################################
  # Calculate mu_j's, Tau_j's for each group:
  muj = matrix(NA, nrow = G, ncol = 2)
  Tauj = array(NA, dim = c(G, 2, 2))
  mj = rep(NA, G)

  for(j in 1:G){

    ###########################################################################
    # Subset jth group data:
    d = df[df$'g' == unique.g[j], ]
    mj[j] = nrow(d)

    ###########################################################################
    # Calculate and match MLE for mean vector:
    muj[j, ] = mle_negbinomial(df = d, eta = MLEeta, lambda = lambda, maxIter = 1000)

    r = muj[j, 1]
    p = muj[j, 2]
    M = mj[j]

    ###########################################################################
    # Calculate and match Hessian values:
    Tauj[j, 1, 2] = Tauj[j, 2, 1] = r * (p - 1) * M
    Tauj[j, 2, 2] = -1 * Tauj[j, 2, 1] + lambda
    Tauj[j, 1, 1] = max(r * M * (r * trigamma(r) + digamma(r) - log(p)) - r * sum(digamma(d$'x' + r) + r * trigamma(d$'x' + r)) + lambda, (0.1 + Tauj[j, 1, 2]^2)/Tauj[j, 2, 2])

  }

  #############################################################################
  # Set muj to appropriate log, log-odds scale, instead of r/p from mle_negbinomial:
  muj[, 2] = muj[, 2] / (1 - muj[, 2])
  muj = log(muj)

  #############################################################################
  # Calculate initial values for mu, Tau; if input 'jitter' is TRUE, then add small jitter to data points:
  if(jitter){
    mu = colSums(jitter(muj) * mj) / sum(mj)
    if(G < 3){
      Tau = solve(crossprod(jitter(muj)))
    } else {
      Tau = solve(cov(jitter(muj)))
    }
  } else {
    mu = colSums(muj * mj) / sum(mj)
    if(G < 3){
      Tau = solve(crossprod(muj))
    } else {
      Tau = solve(cov(muj))
    }
  }

  #############################################################################
  # Initialize empty objects for WHILE loop:
  Grad  = matrix(tol, 2, 2)
  Rhoji = array(NA, dim = c(G, 2, 2))
  Aj    = matrix(NA, nrow = 2, ncol = G)
  iternum = 0

  #############################################################################
  # Calculate initial objective function value:
  B = 0
  for(j in 1:G){

    Rhoji[j, , ] = solve(Tau + Tauj[j, , ])
    C = Tau %*% mu + Tauj[j, , ] %*% muj[j, ]
    B = B + t(C) %*% Rhoji[j, , ] %*% C

  }
  obj = 0.5 * (G * log(det(Tau)) - G * t(mu) %*% Tau %*% mu - sum(log(Rhoji[, 1, 1] * Rhoji[, 2, 2] - Rhoji[, 1, 2] * Rhoji[, 2, 1])) + B)
  err = tol + 1

  #############################################################################
  # Implement while loop that fits Tau/mu by coordinate gradient/analytic descent:
  while((abs(err) > tol) & iternum < maxIter){

    ###########################################################################
    # Prepare calculations for new Tau gradient calculation:
    for(j in 1:G){

      Rhoji[j, , ] = solve(Tau + Tauj[j, , ])
      Aj[, j] = Rhoji[j, , ] %*% (Tau %*% mu + Tauj[j, , ] %*% muj[j, ])

    }

    ###########################################################################
    # Calculate Tau gradient:
    Grad = 0.5 * (G * solve(Tau) - (colSums(Rhoji) + tcrossprod(Aj - mu)))

    ###########################################################################
    # Take Tau step:
    Tau = Tau + EMPBeta * Grad

    ###########################################################################
    # Calculate new mu:
    B = numeric(2)
    for(j in 1:G){

      Rhoji[j, , ] = solve(Tau + Tauj[j, , ])
      B = B + Rhoji[j, , ] %*% Tauj[j, , ] %*% muj[j, ]

    }
    mu = c(solve(diag(2) - colSums(Rhoji) %*% Tau / G) %*% B) / G

    ###########################################################################
    # Update iteration count to exit loop at maxIter:
    iternum = iternum + 1

    ###########################################################################
    # Calculate current objective function value:
    B = 0
    for(j in 1:G){

      C = Tau %*% mu + Tauj[j, , ] %*% muj[j, ]
      B = B + t(C) %*% Rhoji[j, , ] %*% C

    }
    obj_old = obj
    obj = 0.5 * (G * log(det(Tau)) - G * t(mu) %*% Tau %*% mu - sum(log(Rhoji[, 1, 1] * Rhoji[, 2, 2] - Rhoji[, 1, 2] * Rhoji[, 2, 1])) + B)

    ###########################################################################
    # Calculate change in objective function from prior step to current:
    err = obj - obj_old

  }

  #############################################################################
  # Return empirical Bayes solutions, mu and Sigma, and last objective function value for user to modify 'tol' for better convergence:
  return(list('mu' = mu, 'Sigma'  = solve(Tau), 'Obj' = obj))

}

