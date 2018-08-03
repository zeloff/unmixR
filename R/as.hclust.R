#' Hierarchical clustering of classification results
#'
#' Performs a hierarchical clustering of individuals from the distances between
#' pairs.
#' 
#' @param x a matrix with three columns, in which the first two are integers
#' identifying the individuals in each pair and the third is the distance
#' between them.  Such a matrix can easily be obtained by \code{cbind}ing the
#' two objects output by \code{elink()}.
#' @param N the number of individuals
#' @param labels a vector of labels to be assigned to each individual.  Defaults
#' to 1:N
#'
#' @return a \code{hclust} object
#'
#' @seealso \link{elink}
#'
#' @export
as.hclust <- function(x, N, labels = 1:N) {
	# Leaves are negative numbers
	x[x[, 1] <= N, 1] <- -x[x[, 1] <= N, 1]
	x[x[, 2] <= N, 2] <- -x[x[, 2] <= N, 2]

	# Merged clusters are positive and start at 1 (not N + 1)
	x[x[, 1] > N, 1] <- x[x[, 1] > N, 1] - N
	x[x[, 2] > N, 2] <- x[x[, 2] > N, 2] - N

	nz <- nrow(x)
	x <- rbind(x, c(nz, nz - 1, 1))

	dend <- list()
	dend$merge <- as.matrix(x[, 1:2])
	dend$height <- x[, 3]

	ordr <- integer(N)
	next_leaf <- 1

	ff2 <- function(merge_id) {
		pair <- sort(dend$merge[merge_id, ])
		if (all(pair < 0)) {
			ordr[next_leaf + 0:1] <<- -pair
			next_leaf <<- next_leaf + 2
			return()
		} else {
			if (pair[1] < 0) {
				ordr[next_leaf] <<- -pair[1]
				next_leaf <<- next_leaf + 1
				ff2(pair[2])
			} else {
				ff2(pair[1])
				ff2(pair[2])
			}
		}
	}
	
	ff3 <- function(merge_id) {
		pair <- dend$merge[merge_id, ]
		if (pair[1] < 0) {
			ordr[next_leaf] <<- -pair[1]
			next_leaf <<- next_leaf + 1
		} else {
			ff3(pair[1])
		}
		if (pair[2] < 0) {
			ordr[next_leaf] <<- -pair[2]
			next_leaf <<- next_leaf + 1
		} else {
			ff3(pair[2])
		}
	}
	
	ff2(nrow(dend$merge) - 1)
	of2 <- ordr
	ordr <- integer(N)
	next_leaf <- 1
	ff3(nrow(dend$merge) - 1)
	of3 <- ordr

	dend$order <- of3
	dend$labels <- labels
	class(dend) <- "hclust"
	dend
}
