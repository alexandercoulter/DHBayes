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
  Sx = mj = rep(NA, G)
  x = df$'x'
  for(j in 1:G){
    d = x[df$'g' == unique.g[j]]
    Sx[j] = sum(d)
    mj[j] = length(d)
  }
  Sm = sum(mj)

  # Exit if number of groups G is 3 or fewer:
  if(G < 4) stop('Algorithm cannot provide solution for three or fewer groups.')

  # Initialize empty objects for WHILE loop:
  Score  = rep(NA, 3)
  Hessian = matrix(NA, 3, 3)
  iternum = 0

  # Give initial values:
  r = a = b = 1
  rab = c(r, a, b)

  # Calculate a', b':
  a. = a + mj * r
  b. = b + Sx

  # Calculate initial objective function value:
  obj = sum(lgamma(x + r)) + sum(lbeta(a., b.)) - G * lbeta(a, b) - Sm * lgamma(r)

  # Implement while loop that fits Tau/mu by coordinate gradient/analytic descent:
  while((abs(err) > tol) & iternum < maxIter){

    r = rab[1]
    a = rab[2]
    b = rab[3]

    # Calculate a', b':
    a. = a + mj * r
    b. = b + Sx

    # Calculate preparatory values:
    sum_a_digamma = sum(digamma(a.) - digamma(a. + b.))
    sum_b_digamma = sum(digamma(b.) - digamma(a. + b.))
    sum_ma_digamma = sum(mj * (digamma(a.) - digamma(a. + b.)))
    sum_mb_digamma = sum(mj * (digamma(b.) - digamma(a. + b.)))
    sum_x_digamma = sum(digamma(x + r))
    digamma_gaab = G * (digamma(a) - digamma(a + b))
    digamma_gbab = G * (digamma(b) - digamma(a + b))
    sum_a_trigamma = sum(trigamma(a.) - trigamma(a. + b.))
    sum_b_trigamma = sum(trigamma(b.) - trigamma(a. + b.))
    sum_ma_trigamma = sum(mj^2 * (trigamma(a.) - trigamma(a. + b.)))
    sum_mb_trigamma = sum(mj^2 * (trigamma(b.) - trigamma(a. + b.)))
    sum_x_trigamma = sum(trigamma(x + r))
    trigamma_gaab = G * (trigamma(a) - trigamma(a + b))
    trigamma_gbab = G * (trigamma(b) - trigamma(a + b))

    # Calculate Score vector:
    Score = c(sum_ma_digamma + sum_x_digamma - Sm * digamma(r),
              sum_a_digamma - digamma_gaab,
              sum_b_digamma - digamma_gbab)

    # Calculate Hessian components:
    H11 = sum_ma_trigamma + sum_x_trigamma - Sm * trigamma(r)
    H22 = sum_a_trigamma - trigamma_gaab
    H33 = sum_b_trigamma - trigamma_gbab
    H12 = sum_ma_trigamma
    H13 = sum(mj * trigamma(a. + b.))
    H23 = G * trigamma(a + b) - sum(trigamma(a. + b.))

    # Combine Hessian components into Hessian matrix:
    Hessian = matrix(c(H11, H12, H13, H12, H22, H23, H13, H23, H33), 3, 3)

    # Take step:
    Step = solve(Hessian, Score)
    rab = rab - eta * Step

    # Update iteration count to exit loop at maxIter:
    iternum = iternum + 1

    # Calculate current objective function value:
    r = rab[1]
    a = rab[2]
    b = rab[3]
    obj_old = obj
    obj = obj = sum(lgamma(x + r)) + sum(lbeta(a., b.)) - G * lbeta(a, b) - Sm * lgamma(r)

    # Calculate change in objective function from prior step to current:
    err = obj - obj_old

  }

  return(list('r' = r, 'a' = a, 'b' = b))
}

