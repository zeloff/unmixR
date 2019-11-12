#' Prior parameters for gamma function
#'
#' Determines the optimal parameters for the gamma distribution, given the
#' specified sample size and prior distribution, by minimizing the
#' Kullback-Leibler divergence between both.
#'
#' @note Sample size has an upper bound on 170.  The Stirling numbers of the
#' first kind for values greater than that become too large, get translated as
#' infinite, and the optimization function fails.
#'
#' @param n Sample size
#' @param distr The probability distribution for the prior
#' @param ... other parameters passed to the distribution function
#'
#' @return A names vector with parameters \code{a} and \code{b}
#' 
#' @seealso \link{KLD}, \link{stirling1}
#' @examples
#' gamma_prior(50, 'unif')
#' gamma_prior(50, 'pois', lambda = 5)
#' gamma_prior(50, 'lnorm', meanlog = 1.6, sdlog = 0.69)
#'
#' @export
gamma_prior <-
	function(n, distr = c('unif', 'pois', 'nbinom', 'lnorm', 'norm'), ...) {
	
	distr <- match.arg(distr)
	
		if (distr != "unif") {
			g <- do.call(paste0('d', distr), args = list(1:n, ...))
			g <- g/sum(g)
		} else {
			g <- rep(1/n, n)
		}

		stlg <- abs(stirling1(n, n))

		res <- stats::optim(par = c(a = 0.5, b = 0.5), KLD,
												params = list(n = n, g = g, stirling = stlg, d = distr),
												lower = c(1e-3, 1e-4), 
												control = list(maxit = 1e4),
												method = 'L-BFGS-B')

		return(res$par)
	}

