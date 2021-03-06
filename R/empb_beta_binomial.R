#' empb_beta_binomial
#'
#' @param df data.frame object, containing at least columns named 'x' containing non-negative integer values (number of successes), 'n' containing non-negative integer values (number of trials), and 'g' containing group labels.
#' @param eta positive numeric dampening parameter for Newton's method, gradient descent algorithm.
#' @param tol non-negative numeric tolerance parameter for exiting optimization algorithm.
#' @param maxIter positive integer setting maximum number of iterations for optimization algorithm.
#' @param method string controlling optimization method; default 'newton'.
#'
#' @return list object containing empirical Bayes (EMPB) estimates of a, b hyperparameters, assuming df$x ~ binom(p_g, df$n), and p_g ~ beta(a, b), where 'p_g' denotes a group-level parameter.
#' @export
#'
#' @examples
#' # Generate example data:
#' set.seed(31)
#' a = 3
#' b = 9
#'
#' # Number of groups:
#' NG = 10
#'
#' # Creating group IDs:
#' g = replicate(NG, paste(sample(LETTERS, 10), sep="", collapse=""))
#'
#' # Generating 'true' p parameters:
#' p = rbeta(length(g), a, b)
#'
#' # Number of experiments, i.e. rows in df:
#' numexps = 100
#'
#' # Filling df with pseudo data; note the requisite columns 'n', 'x', and 'g':
#' df = data.frame('n' = numeric(0), 'x' = numeric(0), 'g' = character(0))
#' for(k in 1:numexps){
#'   gk = sample(g, 1)
#'   nk = 5 + rpois(1, 10)
#'   xk = rbinom(1, nk, p[g == gk])
#'   df = rbind(df, data.frame('n' = nk, 'x' = xk, 'g' = gk))
#' }
#'
#' # Generating empirical Bayes (EMPB) solutions for a, b:
#' ab_fit = empb_beta_binomial(df = df)
#'
#' # Compare fitted values to known values:
#' cbind(c(a, b), c(ab_fit$a, ab_fit$b))
empb_beta_binomial = function(df, eta = 0.1, tol = 1e-5, maxIter = 10000, method = c('newton', 'gdescent')){

  #############################################################################
  # Object 'df' should be 'data.frame', with columns 'n', 'x', and 'g'.  To that end:
  if(class(df) != 'data.frame') stop('Object \'df\' should be of type \'data.frame\'.')
  if(!all(c('n', 'x', 'g') %in% names(df))) stop('Object \'df\' must contain data.frame columns named \'n\', \'x\', and \'g\'.')
  if((length(df$n) != length(df$x)) | (length(df$n) != length(df$g))) stop('Object \'df\' data.frame columns are not of same length.')

  #############################################################################
  # Calculate parameters, starting points for optimizer:
  unique.g = sort(unique(df$'g'))
  Lg = length(unique.g)
  Ng = Sg = Tg = rep(0, length(unique.g))
  for(j in 1:Lg){

    Sg[j] = sum(df$x[df$g == unique.g[j]])
    Ng[j] = sum(df$n[df$g == unique.g[j]])
    Tg[j] = Ng[j] - Sg[j]

  }

  #############################################################################
  # Exit if any groups have zero sample size:
  if(any(Ng == 0)) stop('Some groups have zero sample size.')

  #############################################################################
  # Initialize values:
  STg = cbind(Sg, Tg)
  ab = c(1, 1)

  #############################################################################
  # If data (df$x) are all zeros, then return(0, mean(Ng)):
  if(sum(Sg) == 0) return(c(0, mean(Ng)))
  # If data (df$x) are all ones, then return(mean(Ng), 0):
  if(sum(Tg) == 0) return(c(mean(Ng), 0))

  #############################################################################
  # Get method:
  method = match.arg(method)

  # Initialize parameter vectors for updating in loop:
  Score = rep(NA, 2)
  Hessian = matrix(NA, 2, 2)
  Step = c(tol, tol)
  iternum = 0

  #############################################################################
  # Run hard-capped Newton's method algorithm:
  while((sum(abs(Step)) > tol) & (iternum < maxIter)){

    ###########################################################################
    # Calculate Score vector:
    Score = Lg * (digamma(sum(ab)) - digamma(ab)) + colSums(digamma(matrix(ab, ncol = 2, nrow = Lg, byrow = TRUE) + STg) - digamma(sum(ab) + Ng))

    ###########################################################################
    # If/else on method input: if 'newton', perform Newton's method; if 'gdescent', perform gradient descent:
    if(method == 'newton'){

      #########################################################################
      # Calculate Hessian matrix:
      Hessian[ , ] = Lg * trigamma(sum(ab)) - sum(trigamma(sum(ab) + Ng))
      Hessian[c(1, 4)] = Hessian[c(1, 4)] - Lg * trigamma(ab) + colSums(trigamma(matrix(ab, ncol = 2, nrow = Lg, byrow = TRUE) + STg))

      #########################################################################
      # Set step:
      Step = solve(Hessian, Score)

    } else {

      #########################################################################
      # Set step:
      Step = -Score

    }

    ###########################################################################
    # Take damped, capped step:
    Step = pmin(Step, 0.9 * ab)
    ab = ab - eta * Step

    ###########################################################################
    # Update iteration count to exit loop at maxIter:
    iternum = iternum + 1

  }

  #############################################################################
  # Return empirical Bayes a, b:
  return(list('a' = ab[1], 'b' = ab[2]))

}
