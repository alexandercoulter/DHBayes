#' mape_poisson
#'
#' @param df data.frame object, containing at least column named 'x' containing non-negative integer values.
#' @param a_prior positive numeric, giving prior parameter 'a', assuming Poisson parameter L ~ beta(a, b).
#' @param b_prior positive numeric, giving prior parameter 'b', assuming Poisson parameter L ~ beta(a, b).
#'
#' @return numeric, maximum a-posteriori estimate (MAPE) of Poisson parameter L, assuming df$x ~ poisson(L), and L ~ beta(a, b).
#' @export
#'
#' @examples
#' # Generate example data:
#' set.seed(31)
#' L = 5
#'
#' # Number of experiments, i.e. rows in df:
#' numexps = 10
#'
#' # Filling df with pseudo data; note the requisite column 'x':
#' df = data.frame('x' = rpois(numexps, L))
#'
#' # Generating maximum a posteriori estimate (MAPE) solution for L:
#' L_fit = mape_poisson(df = df, a_prior = 1, b_prior = 1)
#'
#' # Compare fitted values to known values:
#' cbind(L, L_fit)
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
