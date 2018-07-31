#include "unmixR.h"
// [[Rcpp::depends(RcppArmadillo)]]

// CHECKED

using namespace Rcpp;
using namespace arma;

vec p_for_1(STUD *consts, MNIW *priors, int N, mat data)
{
	vec p(N);
	mat Sigma(priors->lambda_0 * (priors->k_0 + 1.0) / (priors->k_0 * (priors->v_0 - consts->D + 1.0)));
	double v = priors->v_0 - consts->D;
	rowvec mu(priors->mu_0);

	mat inv_Sigma = chol(Sigma);
	double log_det_Sigma = sum(log(diagvec(inv_Sigma)));

	double vd = v + consts->D;
	double d2 = consts->D / 2.0;

	double lpc = consts->pc_gammaln_by2[vd] - (consts->pc_gammaln_by2[v] + d2 * consts->pc_log[v] + d2 * consts->pc_log_pi) - log_det_Sigma;

	for (int i = 0; i < N; i++) {
		rowvec y(data.row(i));
		rowvec u(y - mu);
		vec z = solve(Sigma, trans(u));
		
		double lp = as_scalar(lpc - ((vd + 1.0) / 2.0) * log(1.0 + (1.0 / (v + 1.0)) * dot(u, z)));
		p(i) = lp;
	}
	return p;
}
