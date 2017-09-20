#' Convert exact linkage results to \code{hclust} objects
#'
#' Converts the output of the \code{elink} function to a \code{hclust} object.
#'
#' @seealso \code{\link{elink}}
#' @param x a matrix with each pair of leaves/nodes (first two columns) and
#' the distance between them (third column)
#' @param labels a unique label that identifies each leaf
#'
#' @return A \code{hclust} object
#' @export
as.hclust <- function(x, labels = 1:(nrow(x) + 1)) {
  N <- nrow(x) + 1

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
  print(of2)
  print(of3)

  dend$order <- of3
  dend$labels <- labels
  class(dend) <- "hclust"
  dend
}