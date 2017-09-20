#include "unmixR.h"

using namespace Rcpp;
using namespace arma;

mat rwish(double v, mat S) {

	// Calculate Wishart distribution residuals
	int n = S.n_cols;
	mat x = randn<mat>(n, n);
	for (int i = 0; i < n; i++) {
		x.diag()[i] = sqrt(as<double>(rchisq(1, v - i)));
	}
	x = trimatu(x) * chol(S);
	x = trans(x) * x;

	return x;
}

mat rwishart(double v, mat s) {
	return rwish(v, s);
}

// [[Rcpp::export]]
SEXP rwishart(SEXP v, NumericMatrix s) {
	double V = (as<NumericVector>(v))[0];
	mat S(s.begin(), s.nrow(), s.ncol(), false);
	return wrap(rwish(V, S));
}


mat rinvwishart(double v, mat s) {
	return inv(rwish(v, inv(s)));
}

// [[Rcpp::export]]
SEXP rinvwishart(NumericVector v, NumericMatrix s) {
	double V = (as<NumericVector>(v))[0];
	mat S(s.begin(), s.nrow(), s.ncol(), false);
	return wrap(inv(rwish(V, inv(S))));
}
