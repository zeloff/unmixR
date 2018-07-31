#' Kullback-Leibler divergence between priors
#'
#' Calculates the Kullback-Leibler divergence between a probability mass
#' function \eqn{g(k)} and a \eqn{Gamma(a, b)} prior for the alpha parameter.
#'
#' This function is not meant to be used by itself, but instead to be called by
#' the \code{gamma_prior} function, which tries to find the values for
#' \eqn{a} and \eqn{b} that minimize the Kullback-Leibler divergence for a
#' given \eqn{k} and \eqn{n}.
#'
#' Besides \code{k} and \code{n}, a matrix of unsigned Stirling numbers of the
#' first kind for \eqn{n} elements and \eqn{n} cycles must also be passed as a
#' fixed parameter.
#' 
#' @param ab A vector with the initial values for the Gamma(a, b) prior
#' parameters
#' @param params A list of fixed parameters. See details.
#'
#' @return An atomic numeric vector with the value of the Kullback-Leibler
#' divergence
#' 
#' @references
#' Dorazio, R.M. (2009). On selecting a prior for the precision parameter of
#' Dirichelet process mixture models. \emph{Journal of Statistical Planning and
#' Inference}, \strong{139}, 3384-3390. doi:10.1016/j.jspi.2009.03.009
#' 
#' @seealso \link{gamma_prior}, \link{stirling1}
#'
KLD <- function(ab, params) {
	
	a <- ab[1]
	b <- ab[2]

	# Probability mass function for the prior induced on K
	PMF_int <- function(alpha, pars) {
		# Expression from Dorazio2009 is
		#
		#   alpha^(pars$k + pars$a - 1) * exp(-pars$b * alpha) * gamma(alpha) /
		#   gamma(alpha + pars$n)
		#
		# but gamma() has a somewhat short range that makes it innapropriate
		# for integration up to Inf, so we exp(log()) the entire thing to use
		# lgamma() instead:
		exp((pars$k + pars$a - 1) * log(alpha) - pars$b * alpha + lgamma(alpha) -
				lgamma(alpha + pars$n))
	}

	pi_k <- vector(, params$n)
	for (k in 1:params$n) {
		nint <- stats::integrate(PMF_int, lower = 0, upper = Inf,
														 pars = list(a = a, b = b, k = k, n = params$n),
														 stop.on.error = FALSE, abs.tol = 0.)$value
		s <- params$stirling[params$n + 1, k + 1]

		if(nint == 0 | is.infinite(s)) {
			pi_k[k] <- NA
			next
		}

		pi_k[k] <- nint * s / gamma(a) * b^a
	}

	KL <- sum(params$g * log(params$g / pi_k))

	return(KL)
}

