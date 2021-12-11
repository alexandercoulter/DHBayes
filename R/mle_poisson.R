#' mle_poisson
#'
#' @param df data.frame object, containing at least column named 'x' containing non-negative integer values.
#'
#' @return MLE of Poisson distribution parameter L, assuming df$x ~ poisson(L).
#' @export
#'
#' @examples
#' # Hello world
#' x = 2
mle_poisson = function(df){

  #############################################################################
  # Object 'df' should be 'data.frame', with column 'x'.  To that end:
  if(class(df) != 'data.frame') stop('Object \'df\' should be of type \'data.frame\'.')
  if(!('x' %in% names(df))) stop('Object \'df\' must contain data.frame column named \'x\'.')

  #############################################################################
  # Calculate MLE from input:
  return(mean(df$'x'))

}
