#' mle_poisson
#'
#' @param df
#'
#' @return
#' @export
#'
#' @examples
mle_poisson = function(df){
  # Object 'df' should be 'data.frame', with column 'x'.  To that end:

  if(class(df) != 'data.frame') stop('Object \'df\' should be of type \'data.frame\'.')
  if(!('x' %in% names(df))) stop('Object \'df\' must contain data.frame column named \'x\'.')

  # Calculate MLE from input:
  return(mean(df$'x'))
}
