mape_binomial = function(df, a, b){
  # Object 'df' should be 'data.frame' or 'list' type, with elements 'n' and 'x'.  To that end:

  if(typeof(df) != list) stop('Object \'df\' should be of type \'data.frame\' or \'list\'.')
  if(!all(c('n', 'x') %in% names(df))) stop('Object \'df\' must contain list elements (or data.frame columns) named \'n\' and \'x\', respectively.')
  if(length(df$n) != length(df$x)) stop('List elements \'n\' and \'x\' are not of same length.')

  # Inputs 'a' and 'b' should be non-negative:
  if(a < 0) stop('Prior parameter \'a\' should be non-negative.')
  if(b < 0) stop('Prior parameter \'b\' should be non-negative.')

  # Calculate MAPE from input:

}
