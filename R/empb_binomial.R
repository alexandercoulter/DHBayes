empb_beta_binomial = function(df, groupIDs){
  # Object 'df' should be 'data.frame' or 'list' type, with elements 'n', 'x', and 'g'.  To that end:

  if(typeof(df) != 'list') stop('Object \'df\' should be of type \'data.frame\' or \'list\'.')
  if(!all(c('n', 'x', 'g') %in% names(df))) stop('Object \'df\' must contain list elements (or data.frame columns) named \'n\' and \'x\', respectively.')
  if((length(df$n) != length(df$x)) | (length(df$n) != length(df$g))) stop('List elements are not of same length.')

  # Object 'groupIDs' should be a character vector which contains at least as many unique elements as in 'df$g'.  So:
  ## Check that 'groupIDs' is character type.
  if(typeof(groupIDs) != 'character') stop('Object \'groupIDs\' should be character vector.')
  ## Check  that each element of 'df$g' is contained in 'groupIDs':
  if(!all(df$g %in% groupIDs)) stop('Data object \'df\' contains a group ID not contained in \'groupIDs\' vector.')

  ####################################################################
  # Calculate parameters, starting points for optimizer:
  ug = unique(groupIDs)
  Lg = length(ug)
  Ng = Sg = Tg = rep(0, length(ug))
  for(j in 1:Lg){
    Sg[j] = sum(df$x[df$g == ug[j]])
    Ng[j] = sum(df$n[df$g == ug[j]])
    Tg[j] = Ng[j] - Sg[j]
  }
  STg = cbind(Sg, Tg)
  ab_0 = c(1, 1)
  ab = c(-1, -1)

  # If data (df$x) are all zeros, then return(0, mean(Ng)):
  if(sum(Sg) == 0) return(c(0, mean(Ng)))
  # If data (df$x) are all ones, then return(mean(Ng), 0):
  if(sum(Tg) == 0) return(c(mean(Ng), 0))

  ####################################################################
  # Run hard-capped Newton's method algorithm:

  # Initialize parameter vectors for updating in loop:
  Score = rep(NA, 2)
  Hessian = matrix(NA, 2, 2)

  while(!all(ab_0 == ab)){
    ab_0 = ab

    # Calculate Score vector:
    Score = Lg * (digamma(sum(ab_0)) - digamma(ab_0)) + colSums(digamma(ab_0 + STg) - digamma(sum(ab_0) + Ng))

    # Calculate Hessian vector:
    Hessian[ , ] = g * trigamma(sum(ab_0)) - sum(trigamma(sum(ab_0) + Ng))
    Hessian[c(1, 4)] = Hessian[c(1, 4)] - g * trigamma(ab_0) + colSums(trigamma(ab_0 + STg))

    # Calculate step:
    Step = solve(Hessian, Score)

    # Cap step at proportion of ab_0:
    Step = pmin(Step, 0.75 * ab_0)

    # Update ab:
    ab = ab_0 - Step

  }
  return(ab)
}
