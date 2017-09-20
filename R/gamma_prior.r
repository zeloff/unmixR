#' Prior parameters for gamma function
#'
#' @param n Sample size
#' @param distr The probability distribution for the prior
#' @param ... other parameters passed to the distribution function
#'
#' @return A names vector with parameters \code{a} and \code{b}
#'
#' @export
#'
#' @examples
#' gamma_prior(150, distr = 'pois', lambda = 3)
gamma_prior <-
  function(n, distr = c('unif', 'pois', 'nbinom', 'lnorm', 'norm'), ...) {
  
  distr <- match.arg(distr)
  
  g <- do.call(paste0('d', distr), args = list(1:n, ...))
  g <- g/sum(g)

  stlg <- abs(stirling1(n, n))

  res <- stats::optim(par = c(a = 0.5, b = 0.5), KLD,
                      params = list(n = n, g = g, stirling = stlg),
                      lower = c(1e-3, 1e-4), upper = c(100, 100),
                      method = 'L-BFGS-B')

  return(res$par)
}

