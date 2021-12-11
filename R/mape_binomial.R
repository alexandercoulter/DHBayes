#' mape_binomial
#'
#' @param df data.frame object, containing at least columns named 'x' containing non-negative integer values (number of successes), and 'n' containing non-negative integer values (number of trials).
#' @param a_prior positive numeric, giving prior parameter 'a', assuming binomial parameter p ~ beta(a, b).
#' @param b_prior positive numeric, giving prior parameter 'b', assuming binomial parameter p ~ beta(a, b).
#'
#' @return numeric, maximum a-posteriori estimate (MAPE) of binomial parameter p, assuming df$x ~ binom(p, df$n), and p ~ beta(a, b).
#' @export
#'
#' @examples
#' # Hello world
#' x = 2
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
  denom = a_prior + b_prior + sum(df$n) - 2
  if(denom <= 0) stop('Not enough samples to calculate MAPE using given prior parameters.')

  #############################################################################
  # Calculate MAPE p:
  p_mape = min(max((a_post - 1) / denom, 0), 1)

  #############################################################################
  # Return MAPE p:
  return(p_mape)

}
