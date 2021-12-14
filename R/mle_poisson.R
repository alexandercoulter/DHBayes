#' mle_poisson
#'
#' @param df data.frame object, containing at least column named 'x' containing non-negative integer values.
#'
#' @return MLE of Poisson distribution parameter L, assuming df$x ~ poisson(L).
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
#' # Generating maximum likelihood estimate (MLE) solution for L:
#' L_fit = mle_poisson(df = df)
#'
#' # Compare fitted values to known values:
#' cbind(L, L_fit)
mle_poisson = function(df){

  #############################################################################
  # Object 'df' should be 'data.frame', with column 'x'.  To that end:
  if(class(df) != 'data.frame') stop('Object \'df\' should be of type \'data.frame\'.')
  if(!('x' %in% names(df))) stop('Object \'df\' must contain data.frame column named \'x\'.')

  #############################################################################
  # Calculate MLE from input:
  return(mean(df$'x'))

}
