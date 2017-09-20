#include "unmixR.h"

// CHECKED

using namespace Rcpp;
using namespace arma;

void update_prior(STUD consts, int K_plus, std::vector<NORMslice> stats, ivec counts, MNIW* priors) {

	mat muu = zeros<mat>(consts.D, K_plus + 1);
  vec sums = zeros<vec>(consts.D);
	cube sig = zeros<cube>(consts.D, consts.D, K_plus + 1);
	cube invsig = zeros<cube>(consts.D, consts.D, K_plus + 1);
	mat sumsig = zeros<mat>(consts.D, consts.D);
	
	for(uword k = 0; k <= K_plus; ++k) {

		int n = counts[k];

		// Same code as in getlik
		rowvec mu_n(priors->k_0 / (priors->k_0 + n) * priors->mu_0 + n / (priors->k_0 + n) * stats[k].means);

		double v_n = priors->v_0 + n;
		mat S(stats[k].sum_squares - n * (trans(stats[k].means) * stats[k].means));
		rowvec zm_Y(stats[k].means - priors->mu_0);
		mat lambda_n(priors->lambda_0 + S + priors->k_0 * n / (priors->k_0 + n) * trans(zm_Y) * zm_Y);

		sig.slice(k) = rinvwishart(v_n, lambda_n);

		muu.col(k) = trans(mu_n) + chol(sig.slice(k) / (priors->k_0 + n)) * randn(consts.D);

		invsig.slice(k) = inv(sig.slice(k));
		sums += invsig.slice(k) * muu.col(k);
		sumsig += invsig.slice(k);
	}

	mat meansig = inv(sumsig);
	rowvec mu_0 = trans(trans(solve(sums, trans(meansig))) + chol(meansig / priors->k_0) * randn(consts.D));

	double sum2 = 0.0;
	for (uword k = 0; k <= K_plus; k++) {
		vec mmm(trans(muu.col(k) - trans(mu_0)) * invsig.slice(k) * (muu.col(k) - trans(mu_0)));
		sum2 += mmm(0);
	}

	priors->k_0 = rgamma(1, (K_plus + 1.0 + priors->ak_0) / 2.0)[0] * (sum2 + priors->bk_0) / 2.0;
	priors->mu_0 = mu_0;
}

