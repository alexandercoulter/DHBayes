#' posterior_gamma_poisson
#'
#' @param df data.frame object, containing at least columns named 'x' containing non-negative integer values (number of successes) and 'g' containing group labels.
#' @param ab_prior
#' @param dist
#' @param Nsamp
#' @param pred_n
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
#' # Hello world
#' x = 2
posterior_gamma_poisson = function(df, ab_prior = NULL, dist = c('post', 'pred'), pred_n = NULL, Nsamp = 1000, ...){

  #############################################################################
  # Object 'df' should be 'data.frame', with columns 'x' and 'g'.  To that end:
  if(class(df) != 'data.frame') stop('Object \'df\' should be of type \'data.frame\'.')
  if(!all(c('x', 'g') %in% names(df))) stop('Object \'df\' must contain data.frame columns named \'x\' and \'g\'.')
  if(length(df$x) != length(df$g)) stop('Object \'df\' data.frame columns are not of same length.')
  unique.g = sort(unique(df$'g'))

  #############################################################################
  # 'pred_n' controls the fixed 'n' parameters in the case the user wants to sample from the posterior-predictive distribution.
  dist = match.arg(dist)
  if(dist == 'pred'){

    # if 'pred_n' is NULL, throw an error:
    if(is.null(pred_n)){
      stop('If specifying \'dist\' = \'pred\', must specify either a single \'pred_n\' value (non-negative integer), or vector of such values with length equal to number of unique IDs in \'df$g\'.')
    }
    # if 'pred_n' is NOT length(1) nor length(unique.g), throw an error:
    if(!(length(pred_n) %in% c(1, length(unique.g)))){
      stop('If specifying \'pred_n\', must be either a single value (non-negative integer) or a vector of such values equal to number of unique IDs in \'df$g\'.')
    }
    # If 'pred_n' is specified and is not all non-negative integers, then only if 'dist' is 'pred', throw an error:
    if(any(pred_n %% 1 != 0) | (any(pred_n < 0))){
      stop('If obtaining posterior-predictive distribution samples, must specify non-negative integer values for \'pred_n\'.')
    }

  }

  #############################################################################
  # If 'ab_prior' is NULL, then calculate prior a, b parameters from empirical Bayes; otherwise, extract them:
  if(is.null(ab_prior)){
    ab = empb_gamma_poisson(df = df, ...)
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

  for(j in 1:G){
    Sx[j] = sum(df$x[df$g == unique.g[j]])
    Mg[j] = sum(df$g == unique.g[j])
  }

  #############################################################################
  # Calculate posterior parameter values:
  a. = a + Sx
  b. = b + Mg

  #############################################################################
  # If 'dist' is 'post' (for 'posterior distribution'), then generate samples from posterior distribution; otherwise, if 'dist' is 'pred' (for 'posterior predictive distribution'), then generate samples from that:
  if(dist == 'post'){

    ###########################################################################
    # Generate samples:
    X = NULL

  } else {

    ###########################################################################
    # Generate samples:
    X = NULL

  }

  #############################################################################
  # Set X to a matrix and set its names to group IDs:
  X = matrix(X, nrow = Nsamp, byrow = TRUE)
  colnames(X) = unique.g

  #############################################################################
  # Return samples:
  return(X)

}
