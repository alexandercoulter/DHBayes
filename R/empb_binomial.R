empb_binomial = function(df){
  # Object 'df' should be 'data.frame' or 'list' type, with elements 'n' and 'x'.  To that end:

  if(typeof(df) != 'list') stop('Object \'df\' should be of type \'data.frame\' or \'list\'.')
  if(!all(c('n', 'x') %in% names(df))) stop('Object \'df\' must contain list elements (or data.frame columns) named \'n\' and \'x\', respectively.')
  if(length(df$n) != length(df$x)) stop('List elements \'n\' and \'x\' are not of same length.')

  ####################################################################
  # Calculate parameters, starting points for optimizer:
  P = df$x / df$n
  S = sum(df$x)
  N = sum(df$n)
  ST = c(S, N - S)

  # If data (df$x) are all zeros, then return(0, N):
   if(S == 0) return(list('a' = 0, 'b' = N))
  # If data (df$x) are all ones, then return(N, 0):
   if(S == 0) return(list('a' = N, 'b' = 0))
  # Otherwise, initalize with method of moments estimators as if group-level P ~ beta(a_start, b_start):
  mP = mean(P)
  vP = var(P)
  a_start = mP * (mP * (1 - mP) / vP - 1)
  b_start = a_start / mP * (1 - mP)

  ####################################################################
  # Run hard-capped Newton's method algorithm:

  # Initialize parameter vectors for updating in loop:
  ab_0 = c(-1, -1)
  ab = c(a_start, b_start)

  # Initialize score vector, Hessian matrix:
  Score = numeric(2)
  Hessian = matrix(NA, 2, 2)

  # Implement loop where exit condition is equality:
  while(!all(ab_0 == ab)){
    ab_0 = ab

    # Calculate score vector components:
    Score = digamma(ab_0 + ST) + digamma(sum(ab_0)) - digamma(sum(ab_0) + N) - digamma(ab_0)

    # Calculate Hessian matrix components:
    Hessian[1:2, 1:2] = trigamma(sum(ab_0)) - trigamma(sum(ab_0) + N)
    Hessian[c(1, 4)] = Hessian[c(1, 4)] + trigamma(ab_0 + ST) - trigamma(ab_0)

    # Calculate step:
    Step = solve(Hessian, Score)

    # Cap step size to some fraction of distance from edge, to prevent over-shooting:
    # Taking fraction to be 0.75:
    Step = pmin(Step, 0.75 * ab_0)

    # Take the step:
    ab = ab_0 - Step
  }
  return(ab)
}
