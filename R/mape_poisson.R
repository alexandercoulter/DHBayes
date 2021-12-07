#' mape_poisson
#'
#' @param df
#' @param a_prior
#' @param b_prior
#'
#' @return
#' @export
#'
#' @examples
mape_poisson = function(df, a_prior, b_prior){

  #############################################################################
  # Object 'df' should be 'data.frame', with column 'x'.  To that end:
  if(class(df) != 'data.frame') stop('Object \'df\' should be of type \'data.frame\'.')
  if(!('x' %in% names(df))) stop('Object \'df\' must contain data.frame column named \'x\'.')

  #############################################################################
  # Compatibility checks:
  # Prior parameters a_prior, b_prior must be positive:
  if(a_prior <= 0) stop('Object \'a_prior\' must be positive.')
  # Prior parameter p_prior must be strictly between 0 and 1:
  if(b_prior <= 0) stop('Object \'b_prior\' must be positive.')

  #############################################################################
  # Calculate a_post, b_post:
  a_post = a_prior + sum(df$'x')
  b_post = b_prior + length(df$'x')

  #############################################################################
  # Return MAPE L(ambda):
  L = max(0, (a_post - 1) / b_post)
  return(L)

}
