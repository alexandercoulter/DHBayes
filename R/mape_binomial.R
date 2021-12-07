#' mape_binomial
#'
#' @param a_prior
#' @param b_prior
#' @param df
#'
#' @return
#' @export
#'
#' @examples
mape_binomial = function(df, a_prior, b_prior){

  #############################################################################
  # Object 'df' should be 'data.frame', with columns 'x' and 'n'.  To that end:
  if(class(df) != 'data.frame') stop('Object \'df\' should be of type \'data.frame\'.')
  if(!all(c('n', 'x') %in% names(df))) stop('Object \'df\' must contain data.frame columns named \'x\' and \'n\'.')
  if(length(df$n) != length(df$x)) stop('Object \'df\' data.frame columns \'n\' and \'x\' are not of same length.')

  #############################################################################
  # Inputs 'a_prior' and 'b_prior' should be non-negative:
  if(a_prior < 0) stop('Prior parameter \'a_prior\' should be non-negative.')
  if(b_prior < 0) stop('Prior parameter \'b_prior\' should be non-negative.')

  #############################################################################
  # Calculate MAPE from input:
  # Calculate numerator:
  a_post = a_prior + sum(df$x)
  # Calculate denominator; if too small, then exit with message:
  denom = a_prior + b + sum(df$n) - 2
  if(denom <= 0) stop('Not enough samples to calculate MAPE using given prior parameters.')

  #############################################################################
  # Calculate MAPE p:
  p_mape = min(max((a_post - 1) / denom, 0), 1)

  #############################################################################
  # Return MAPE p:
  return(p_mape)

}
