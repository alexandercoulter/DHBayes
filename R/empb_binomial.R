empb_beta_binomial = function(df, groupIDs, eta = 0.1, tol = 1e-5){
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
  Lg = sum(as.numeric(Ng > 0))
  STg = cbind(Sg, Tg)
  ab_0 = c(-1, -1)
  ab = c(1, 1)

  # If data (df$x) are all zeros, then return(0, mean(Ng)):
  if(sum(Sg) == 0) return(c(0, mean(Ng)))
  # If data (df$x) are all ones, then return(mean(Ng), 0):
  if(sum(Tg) == 0) return(c(mean(Ng), 0))

  ####################################################################
  # Run hard-capped Newton's method algorithm:

  # Initialize parameter vectors for updating in loop:
  Score = rep(NA, 2)
  Hessian = matrix(NA, 2, 2)
  Step = 10 * c(tol, tol)

  while(sum(abs(Step)) > tol){

    # Calculate Score vector:
    Score = Lg * (digamma(sum(ab)) - digamma(ab)) + colSums(digamma(ab + STg) - digamma(sum(ab) + Ng))

    # Calculate Hessian vector:
    Hessian[ , ] = Lg * trigamma(sum(ab)) - sum(trigamma(sum(ab) + Ng))
    Hessian[c(1, 4)] = Hessian[c(1, 4)] - Lg * trigamma(ab) + colSums(trigamma(ab + STg))

    # Calculate step:
    Step = eta * solve(Hessian, Score)

    # Cap step at proportion of ab_0:
    Step = pmin(Step, 0.9 * ab)

    # Update ab:
    ab = ab - Step
    print(sum(abs(Step)))
  }
  return(ab)
}
