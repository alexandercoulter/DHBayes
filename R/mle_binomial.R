mle_binomial = function(df){
  # Object 'df' should be 'data.frame' or 'list' type, with elements 'n' and 'x'.  To that end:

  df = data.frame('n' = c(10, 10, 10),
                  'x' = c(3, 4, 3))
  dfL = list('n' = c(10, 10, 10),
             'x' = c(3, 4, 4))
  if(typeof(dfL) != list) stop('Object \'df\' should be of type \'data.frame\' or \'list\'.')
  if(!all(c('n', 'x') %in% names(dfL))) stop('Object \'df\' must contain list elements (or data.frame columns) named \'n\' and \'x\', respectively.')
  if(length(df$n) != length(df$x)) stop('List elements \'n\' and \'x\' are not of same length.')

}
