#include "unmixR.h"

// CHECKED

using namespace arma;

void preclass(STUD *consts, std::vector<NORMslice>* stats, MNIW *priors, ivec *allcounts, ivec *counts, ivec *class_id) {

	vec distances(consts->nsources);

	for (int i = 0; i < consts->N; ++i) {

		for(int j = 0; j < consts->nsources; ++j) {
			distances(j) = sqrt(sum(pow(consts->datas.row(i) - trans(consts->orgmeans.col(j)), 2)));
		}	

		// Randomly pick a source 
		vec cs = 1 - cumsum(distances / sum(distances));
		
		double rand = as_scalar(randu(1));
// TODO: REMOVE THIS
//		double rand = 0.5;
		int choose = as_scalar(find(cs < rand, 1));

		(*allcounts)(choose)++;
		(*counts)(choose)++;
		(*class_id)(i) = choose;

		int counts = (*allcounts)(choose);

		// Recalculate that source's stats
		(*stats)[choose].means += 1.0 / counts * (consts->datas.row(i) - (*stats)[choose].means);
		(*stats)[choose].sum_squares += trans(consts->datas.row(i)) * consts->datas.row(i);

		double k_n = priors->k_0 + counts;
		double v_n = priors->v_0 + counts;

		rowvec zm_Y = (*stats)[choose].means - priors->mu_0;

		mat lambda_n = priors->lambda_0 +
			((*stats)[choose].sum_squares - counts * (trans((*stats)[choose].means) * (*stats)[choose].means)) +
			priors->k_0 * counts / k_n *
			(trans(zm_Y) * zm_Y);

		mat Sigma = lambda_n * (k_n + 1) / (k_n * (v_n - consts->D + 1));

		mat inv_cov = chol(Sigma, "upper");
		(*stats)[choose].inv_cov = inv_cov;
		(*stats)[choose].log_det_cov = sum(log(inv_cov.diag()));
	}
}

