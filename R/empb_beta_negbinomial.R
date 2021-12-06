empb_beta_negbinomial = function(df, eta = 0.01, tol = 0.1, maxIter = 10000){
  # Object 'df' should be 'data.frame' or 'list' type, with elements 'x' and 'g'.  To that end:

  if(typeof(df) != 'list') stop('Object \'df\' should be of type \'data.frame\' or \'list\'.')
  if(!all(c('x', 'g') %in% names(df))) stop('Object \'df\' must contain list elements (or data.frame columns) named \'x\' and \'g\', respectively.')
  if(length(df$x) != length(df$g)) stop('List elements are not of same length.')

  # Dampening parameter eta must be positive:
  if(eta <= 0) stop('Object \'eta\' must be positive.')
  # Tolerance must be non-negative:
  if(tol < 0) stop('Object \'tol\' must be non-negative.')
  # Maximum number of iterations must be positive integer:
  if((maxIter < 1) | ((maxIter %% 1) != 0)) stop('Object \'maxIter\' must be positive integer.')

  # Pull out unique group IDs:
  unique.g = sort(unique(df$'g'))
  G = length(unique.g)

  # Exit if number of groups G is 3 or fewer:
  if(G < 4) stop('Algorithm cannot provide solution for three or fewer groups.')

  # Initialize empty objects for WHILE loop:
  Score  = rep(NA, 3)
  Hessian = matrix(NA, 3, 3)
  rab = rep(NA, 3)
  iternum = 0

  # Calculate initial objective function value:
  obj = NULL

  # Implement while loop that fits Tau/mu by coordinate gradient/analytic descent:
  while((abs(err) > tol) & iternum < maxIter){

    # Calculate a', b':


    # Calculate Score vector:


    # Calculate Hessian components:


    # Combine Hessian components into Hessian matrix:


    # Take step:


    # Update iteration count to exit loop at maxIter:
    iternum = iternum + 1

    # Calculate current objective function value:
    obj_old = obj
    obj = NULL

    # Calculate change in objective function from prior step to current:
    err = obj - obj_old

  }

  return(list('r' = rab[1], 'a' = rab[2], 'b' = rab[3]))
}

