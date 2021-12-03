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
empb_norm_negbinomial = function(df, lambda = 0.01, MLEeta = 0.1, EMPBeta = 0.001, tol = 1e-5, maxIter = 200){
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
  unique.g = unique(df$'g')
  G = length(unique.g)

  # Calculate mu_j's for each group:
  muj = matrix(NA, nrow = G, ncol = 2)
  for(j in 1:G) muj[j, ] = mle_negbinomial(df = df[df$'g' == unique.g[j], ], eta = MLEeta, lambda = lambda, tol = tol, maxIter = maxIter)

  return(list('muj' = muj, 'groups' = unique.g))
  # Calculate Tau_j's for each group:
  Tauj = array(NA, dim = c(G, 2, 2))

  # Calculate initial values for mu, Tau:
  mu  = NULL
  Tau = NULL

  # Initialize empty objects for WHILE loop:
  Grad = matrix(NA, 2, 2)
  Rhoj = array(NA, dim = c(G, 2, 2))
  Aj   = array(NA, dim = c(G, 2, 2))
  iternum = 0

  while((sum(abs(Grad)) > tol) & iternum < maxIter){

    # Prepare calculations for new Tau gradient calculation:


    # Calculate Tau gradient:
    Grad = NULL

    # Take Tau step:
    Tau = Tau - eta * Grad

    # Calculate new mu:
    mu = NULL

    iternum = iternum + 1
    if((det(Grad) < 0) | (det(Grad[1, 1, drop = FALSE]) < 0)) stop("Step isn't positive semi-definite!")
  }

  return(list('mu' = mu, 'Sigma'  = solve(Tau)))
}
