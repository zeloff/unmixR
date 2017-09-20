#include "unmixR.h"

// CHECKED

//# define D(X) do {} while (0)

//using namespace Rcpp;
using namespace arma;

void preallocate(STUD *consts, std::vector<NORMslice>* stats, MNIW *priors, ivec *allcounts) {

	for(int i = 0; i < consts->nsources; ++i) {

		int counts;
		int src;
		rowvec means(consts->D);
		mat ss1(consts->D, consts->D);
		mat ss2(consts->D, consts->D);

		src = consts->sources[i];
		uvec idxs = find(consts->label == src);

		counts = idxs.n_elem;

		means = mean(consts->baseline.rows(idxs), 0);

		ss1 = cov(consts->baseline.rows(idxs)) * (counts - 1);
		ss2 = ss1 + counts * (trans(means) * means);

		double k_n = priors->k_0 + counts;
		double v_n = priors->v_0 + counts;

		rowvec zm_Y = means - priors->mu_0;

		mat lambda_n = priors->lambda_0 + ss1 +
			priors->k_0 * counts / (priors->k_0 + counts) *
			(trans(zm_Y) * zm_Y);

		mat Sigma = lambda_n * (k_n + 1) / (k_n * (v_n - consts->D + 1));

		mat inv_cov = chol(Sigma, "upper");
		double log_det_cov = sum(log(inv_cov.diag()));
		/*

	double vd = v + consts->D;
	double d2 = consts->D / 2.0;

		try {
			stats->inv_cov = chol(Sigma, "upper");
			stats->log_det_cov = sum(log(diagvec(stats->inv_cov)));
		} catch ( std::exception& __ex__ ) {
			throw __ex__;
		}
*/

		consts->ns(src) = counts;
		consts->orgmeans.col(src) = trans(means);
		consts->orgsum_squares.slice(src) = ss2;
		consts->orginv_cov.slice(src) = inv_cov;
		consts->orglog_det_cov(src) = log_det_cov;

		(*stats)[src].means = means;
		(*stats)[src].sum_squares = ss2;
		(*stats)[src].inv_cov = inv_cov;
		(*stats)[src].log_det_cov = log_det_cov;

		(*allcounts)(src) = counts;
	}
}

