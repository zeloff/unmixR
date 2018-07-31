#' Prior parameters for gamma function
#'
#' @param n Sample size
#' @param distr The probability distribution for the prior
#' @param ... other parameters passed to the distribution function
#'
#' @return A names vector with parameters \code{a} and \code{b}
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

