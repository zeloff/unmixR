#include "unmixR.h"

// CHECKED

void update_alpha(double* alpha, int N, int K_plus, MNIW priors) {
	double nu = Rcpp::rbeta(1, *alpha + 1.0, N)[0];
	double bln = priors.b_0 - log(nu);
	double pis = (priors.a_0 + K_plus) / (priors.a_0 + K_plus + N * bln);
	double rg1 = Rcpp::rgamma(1, priors.a_0 + K_plus + 1)[0]; 
	double rg2 = Rcpp::rgamma(1, priors.a_0 + K_plus)[0]; 
	*alpha = pis * rg1 / bln + (1.0 - pis) * rg2 / bln;
}
