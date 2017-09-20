#include "unmixR.h"
#include <math.h>

// CHECKED

using namespace Rcpp;
using namespace arma;

double getlik(STUD *consts, MNIW *priors, NORMslice *stats, rowvec y, int n, bool lik, bool suffs) {

	rowvec mu(priors->k_0 / (priors->k_0 + n) * priors->mu_0 + n / (priors->k_0 + n) * stats->means);

	double k_n = priors->k_0 + n;
	double v_n = priors->v_0 + n;

	mat S(stats->sum_squares - n * (trans(stats->means) * stats->means));
	rowvec zm_Y(stats->means - priors->mu_0);
	mat lambda_n(priors->lambda_0 + S + priors->k_0 * n / (priors->k_0 + consts->N) * trans(zm_Y) * zm_Y);

	double v = v_n - consts->D + 1;
	mat Sigma(lambda_n * (k_n + 1) / (k_n * v));

	double vd = v + consts->D;
	double d2 = consts->D / 2.0;

	if (suffs) {
		try {
			stats->inv_cov = chol(Sigma, "upper");
			stats->log_det_cov = sum(log(diagvec(stats->inv_cov)));
		} catch (...) {
			D("[getlik]   stats(means): " << stats->means);
			D("[getlik]    priors(k_0):   " << priors->k_0 << "\n");
			D("[getlik]   priors(mu_0): " << priors->mu_0);
			D("[getlik]              y: " << y);
			D("[getlik]              n:    " << n << "\n");
			D("[getlik]       lambda_n:\n" << lambda_n);
			D("[getlik]          Sigma:\n" << Sigma);
			D("[getlik]        inv_cov:\n" << stats->inv_cov);
			if (lik) {
				throw;
			}
		}
	}

	double lp, lpx;

	if (!(suffs) | lik) {
		rowvec u(y - mu);
		mat z = solve(stats->inv_cov, trans(u));
		lpx = as_scalar(consts->pc_gammaln_by2(vd - 1) - (consts->pc_gammaln_by2(v - 1) + d2 * (consts->pc_log(v - 1) + consts->pc_log_pi)) - stats->log_det_cov - (vd / 2.0) * log(1.0 + (1.0 / v) * dot(u, z)));

		v = priors->v_0 + n;
		vd = v + consts->D;

		double p = as_scalar(tgamma(vd / 2.0) / (tgamma(v / 2.0) * pow(v, consts->D / 2.0) * pow(M_PI, consts->D / 2.0)) * sqrt(exp(stats->log_det_cov)) * pow(1.0 + (1.0 / v) * (u * stats->inv_cov * trans(u)), -(vd / 2.0)));
		lp = log(p);
		
		D("[getlik] lp: " << lp << "  lpx: " << lpx << "\n");

		if (ISNAN(lp)) {
			D("[getlik]              n: " << n << "\n");
			D("[getlik]            v_0: " << priors->v_0 << "\n");
			D("[getlik]            v_n: " << v_n << "\n");
			D("[getlik]              v: " << v << "\n");
			D("[getlik]             vd: " << vd << "\n");
			D("[getlik]              u: " << u);
			D("[getlik]              z: " << trans(z));
			D("[getlik]            dot: " << dot(u, z) << "\n");
			D("[getlik]            log: " << 1.0 + (1.0 / v) * dot(u, z) << "\n");
			throw std::domain_error("Invalid likelihood. Stopping simulation.\n");
		}
		return lp;
	} 
	return -1; // should be ignored up the stack
}
