#' mle_binomial
#'
#' @param df
#'
#' @return
#' @export
#'
#' @examples
mle_binomial = function(df){
  # Object 'df' should be 'data.frame', with columns 'x' and 'n'.  To that end:

  if(class(df) != 'data.frame') stop('Object \'df\' should be of type \'data.frame\'.')
  if(!all(c('n', 'x') %in% names(df))) stop('Object \'df\' must contain data.frame columns named \'x\' and \'n\'.')
  if(length(df$n) != length(df$x)) stop('Object \'df\' data.frame columns \'n\' and \'x\' are not of same length.')

  # Calculate MLE from input:
  return(sum(df$x) / sum(df$n))
}
