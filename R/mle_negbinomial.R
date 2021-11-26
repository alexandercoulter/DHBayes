mle_negbinomial = function(df, eta = 1, tol = 0.0001){
  # Object 'df' should be 'data.frame' or 'list' type, with element 'x'.  To that end:

  if(typeof(df) != 'list') stop('Object \'df\' should be of type \'data.frame\' or \'list\'.')
  if(!('x' %in% names(df))) stop('Object \'df\' must contain list element (or data.frame column) named \'x\'.')

  # Calculate MLE from input:
  M = length(df$x)
  Sx = sum(df$x)
  MSx = M / Sx

  r = 1
  o = MSx
  rstep = tol + 1
  while(abs(rstep) > tol){
    # Calculate derivative for Newton's step:
    f. = M * (log(o / (1 + o)) - digamma(r)) + sum(digamma(df$x + r))

    # Calculate second derivative for Newton's step:
    f.. = sum(trigamma(df$x + r)) - M * trigamma(r)

    # Calculate step for 'r', capping at some proportion of existing 'r' value to not overshoot into negatives:
    rstep = min(0.5 * r, f. / f..)

    # Take damped step:
    r = r - eta * rstep

    # Calculate maximizing 'o':
    o = r * MSx
  }
  return(c(r, o / (1 + o)))
}
