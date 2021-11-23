empb_binomial = function(df){
  # Object 'df' should be 'data.frame' or 'list' type, with elements 'n' and 'x'.  To that end:

  if(typeof(df) != list) stop('Object \'df\' should be of type \'data.frame\' or \'list\'.')
  if(!all(c('n', 'x') %in% names(df))) stop('Object \'df\' must contain list elements (or data.frame columns) named \'n\' and \'x\', respectively.')
  if(length(df$n) != length(df$x)) stop('List elements \'n\' and \'x\' are not of same length.')

  ####################################################################
  # Calculate parameters, starting points for optimizer:
  P = df$x / df$n
  S = sum(df$x)
  N = sum(df$n)

  # If data (df$x) are all zeros, then return(0, N):
   if(S == 0) return(list('a' = 0, 'b' = N))
  # If data (df$x) are all ones, then return(N, 0):
   if(S == 0) return(list('a' = N, 'b' = 0))
  # Otherwise, initalize with method of moments estimators as if group-level P ~ beta(a_start, b_start):
  mP = mean(P)
  vP = var(P)
  a_start = mP * (mP * (1 - mP) / vP - 1)
  b_start = a_start / mP * (1 - mP)

  return(p_mape)
}
