#' posterior_beta_negbinomial
#'
#' @param df data.frame object, containing at least columns named 'x' containing non-negative integer values and 'g' containing group labels.
#' @param rab_prior length-3, positive numeric vector specifying prior hyperparameter values for r, a, and b; if NULL, values fit by empirical Bayes (EMPB).
#' @param dist string specifying what type of samples to return: either 'post' (for samples from the posterior distribution), or 'pred' (for samples from the posterior-predictive distribution).
#' @param Nsamp positive integer, number of samples to generate per group.
#' @param ... optional parameters to be passed to control EMPB convergence, in the case 'rab_prior' is NULL; see 'empb_gamma_poisson'.
#'
#' @return matrix of samples from the posterior (or posterior-predictive) distribution, where (named) columns are for group IDs included in df$g, and rows are samples; assuming df$x ~ nbinom(r, p_g), and p_g ~ beta(a, b), where 'p_g' denotes a group-level parameter.
#' @export
#'
#' @examples
#' # Hello world
#' x = 2
posterior_beta_negbinomial = function(df, rab_prior = NULL, dist = c('post', 'pred'), Nsamp = 1000, ...){

  #############################################################################
  # Object 'df' should be 'data.frame', with columns 'x' and 'g'.  To that end:
  if(class(df) != 'data.frame') stop('Object \'df\' should be of type \'data.frame\'.')
  if(!all(c('x', 'g') %in% names(df))) stop('Object \'df\' must contain data.frame columns named \'x\' and \'g\'.')
  if(length(df$x) != length(df$g)) stop('Object \'df\' data.frame columns are not of same length.')
  unique.g = sort(unique(df$'g'))

  #############################################################################
  # If 'rab_prior' is NULL, then calculate prior r, a, and b parameters from empirical Bayes; otherwise, extract them:
  if(is.null(rab_prior)){
    rab = empb_beta_negbinomial(df = df, ...)
    r = rab$r
    a = rab$a
    b = rab$b
  } else {
    # rab_prior must be length 3, positive:
    if((length(rab_prior) != 3) | any(rab_prior < 0)) stop('If specifying \'rab_prior\', must be a length-3, strictly positive numeric vector.')
    r = rab_prior[1]
    a = rab_prior[2]
    b = rab_prior[3]
  }

  #############################################################################
  # Calculate group-level data for posterior parameters:
  unique.g = sort(unique(df$'g'))
  G = length(unique.g)
  Sx = Mg = rep(0, G)

  for(j in 1:G){
    Sx[j] = sum(df$x[df$g == unique.g[j]])
    Mg[j] = sum(df$g == unique.g[j])
  }

  #############################################################################
  # Calculate posterior parameter values:
  a. = a + r * Mg
  b. = b + Sx

  #############################################################################
  # Get distribution to sample from:
  dist = match.arg(dist)

  #############################################################################
  # If 'dist' is 'post' (for 'posterior distribution'), then generate samples from posterior distribution; otherwise, if 'dist' is 'pred' (for 'posterior predictive distribution'), then generate samples from that:
  X = rbeta(n = Nsamp * G,
            shape1 = a.,
            shape2 = b.)

  if(dist == 'pred'){
    X = rnbinom(n = Nsamp * G,
                size = r,
                prob = X)
  }

  #############################################################################
  # Set X to a matrix and set its names to group IDs:
  X = matrix(X, nrow = Nsamp, byrow = TRUE)
  colnames(X) = unique.g

  #############################################################################
  # Return samples:
  return(X)

}
