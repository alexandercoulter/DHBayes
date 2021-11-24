#' mle_binomial
#'
#' @param df
#'
#' @return
#' @export
#'
#' @examples
mle_binomial = function(df){
  # Object 'df' should be 'data.frame' or 'list' type, with elements 'n' and 'x'.  To that end:

  if(typeof(df) != 'list') stop('Object \'df\' should be of type \'data.frame\' or \'list\'.')
  if(!all(c('n', 'x') %in% names(df))) stop('Object \'df\' must contain list elements (or data.frame columns) named \'n\' and \'x\', respectively.')
  if(length(df$n) != length(df$x)) stop('List elements \'n\' and \'x\' are not of same length.')

  # Calculate MLE from input:
  return(sum(df$x) / sum(df$n))
}
