#' Stirling numbers of the first kind
#' 
#' Calculation of (signed) Stirling numbers of the first kind for a given number
#' of elements and orbits.
#'
#' @param n Number of elements.
#' @param k Number of orbits.
#'
#' @return A lower triangular matrix of Stirling numbers of the first kind for
#' every pair of (n, k).
#'  
#' @export
#' @examples
#' stirling1(4, 5)
#'
stirling1 <- function(n, k) {
	if (!is.numeric(k) | k < 0) stop("k must be a non-negative integer")
	if (!is.numeric(n) | n < 0) stop("n must be a non-negative integer")

	mm <- diag(1, n + 1, k + 1)
	mm[upper.tri(mm)] <- NA

	if (n < 2 | k < 1) return(mm)

	for (nn in seq.int(2, n)) {
		maxk <- pmin(k + 1, nn)
		mm[nn + 1, 2:maxk] <- -(nn - 1) * mm[nn, 2:maxk] + mm[nn, 1:(maxk - 1)]
	}

	dimnames(mm) <- list(n = 0:n, k = 0:k)
	return(mm)	
}

