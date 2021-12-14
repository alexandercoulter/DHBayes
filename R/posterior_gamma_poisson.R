#' posterior_gamma_poisson
#'
#' @param df data.frame object, containing at least columns named 'x' containing non-negative integer values and 'g' containing group labels.
#' @param ab_prior length-2, positive numeric vector specifying prior hyperparameter values for a, b; if NULL, values fit by empirical Bayes (EMPB).
#' @param dist string specifying what type of samples to return: either 'post' (for samples from the posterior distribution), or 'pred' (for samples from the posterior-predictive distribution).
#' @param Nsamp positive integer, number of samples to generate per group.
#' @param ... optional parameters to be passed to control EMPB convergence, in the case 'ab_prior' is NULL; see 'empb_gamma_poisson' or 'empb_gamma_poisson_c'.
#'
#' @return matrix of samples from the posterior (or posterior-predictive) distribution, where (named) columns are for group IDs included in df$g, and rows are samples; assuming df$x ~ poisson(L_g), and L_g ~ gamma(a, b), where 'L_g' denotes a group-level parameter.
#' @export
#'
#' @examples
#' # Generate example data:
#' set.seed(31)
#' a = 23
#' b = 9
#'
#' # Number of groups:
#' NG = 10
#'
#' # Creating group IDs:
#' g = replicate(NG, paste(sample(LETTERS, 10), sep="", collapse=""))
#'
#' # Generating 'true' L parameters:
#' L = rgamma(length(g), a, b)
#'
#' # Number of experiments, i.e. rows in df:
#' numexps = 100
#'
#' # Filling df with pseudo data; note the requisite columns 'x' and 'g':
#' df = data.frame('x' = numeric(0), 'g' = character(0))
#' for(k in 1:numexps){
#'   gk = sample(g, 1)
#'   xk = rpois(1, L[g == gk])
#'   df = rbind(df, data.frame('x' = xk, 'g' = gk))
#' }
#'
#' # Generating 1000 posterior distribution samples for each group:
#' posterior_values = posterior_gamma_poisson(df = df)
#' dim(posterior_values)
#'
#' # Create histogram of posterior distribution samples for first group (by alphabetic order):
#' hist(posterior_values[, 1], main = colnames(posterior_values)[1])
posterior_gamma_poisson = function(df, ab_prior = NULL, dist = c('post', 'pred'), Nsamp = 1000, ...){

  #############################################################################
  # Object 'df' should be 'data.frame', with columns 'x' and 'g'.  To that end:
  if(class(df) != 'data.frame') stop('Object \'df\' should be of type \'data.frame\'.')
  if(!all(c('x', 'g') %in% names(df))) stop('Object \'df\' must contain data.frame columns named \'x\' and \'g\'.')
  if(length(df$x) != length(df$g)) stop('Object \'df\' data.frame columns are not of same length.')

  #############################################################################
  # If 'ab_prior' is NULL, then calculate prior a, b parameters from empirical Bayes; otherwise, extract them:
  if(is.null(ab_prior)){
    ab = empb_gamma_poisson_c(df = df, ...)
    a = ab$a
    b = ab$b
  } else {
    # ab_prior must be length 2, positive:
    if((length(ab_prior) != 2) | any(ab_prior < 0)) stop('If specifying \'ab_prior\', must be a length-2, strictly positive numeric vector.')
    a = ab_prior[1]
    b = ab_prior[2]
  }

  #############################################################################
  # Calculate group-level data for posterior parameters:
  unique.g = sort(unique(df$'g'))
  G = length(unique.g)
  Sx = Mg = rep(0, G)
  x = df$x
  g = df$g

  for(j in 1:G){
    d = x[g == unique.g[j]]
    Sx[j] = sum(d)
    Mg[j] = length(d)
  }

  #############################################################################
  # Calculate posterior parameter values:
  a. = a + Sx
  b. = b + Mg

  #############################################################################
  # Get distribution to sample from:
  dist = match.arg(dist)

  #############################################################################
  # If 'dist' is 'post' (for 'posterior distribution'), then generate samples from posterior distribution; otherwise, if 'dist' is 'pred' (for 'posterior predictive distribution'), then generate samples from that:
  if(dist == 'post'){

    ###########################################################################
    # Generate samples:
    X = rgamma(n = Nsamp * G,
               shape = a.,
               rate = b.)

  } else {

    ###########################################################################
    # Generate samples:
    X = rnbinom(n = Nsamp * G,
                size = a.,
                prob = b. / (1 + b.))

  }

  #############################################################################
  # Set X to a matrix and set its names to group IDs:
  X = matrix(X, nrow = Nsamp, byrow = TRUE)
  colnames(X) = unique.g

  #############################################################################
  # Return samples:
  return(X)

}
