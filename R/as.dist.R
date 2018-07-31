as.dist <- function(x, labels) {
	ff <- function(a, b, mat) {
	
		get_dist <- function(x, y) {
			sum(mat[x, ] != mat[y, ]) / ncol(mat)
		}

		mapply(get_dist, a, b)
	}

	out <- outer(1:nrow(x), 1:nrow(x), ff, x)
	dimnames(out) <- list(labels, labels)
	as.dist(out)
}
