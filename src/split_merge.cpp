#include "unmixR.h"
// # define D(X) do {} while (0)
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

void split_merge(mat data, ivec *source_id, STUD *consts, int N, MNIW *priors, std::vector<NORMslice>* stats, ivec *counts, ivec *allcounts, int *K_plus, double alpha)
{
	// INITIALIZE STUFF
	
	// Individual stats
	std::vector<NORMslice> sstats(2);
	for (int d = 0; d < 2; d++) {
			sstats[d].means = zeros<rowvec>(consts->D);
			sstats[d].sum_squares = zeros<mat>(consts->D, consts->D);
			sstats[d].inv_cov = zeros<mat>(consts->D, consts->D);
			sstats[d].log_det_cov = zeros<vec>(1);
	}

	// Combined stats
	NORMslice cstats;
	cstats.means = zeros<rowvec>(consts->D);
	cstats.sum_squares = zeros<mat>(consts->D, consts->D);
	cstats.inv_cov = zeros<mat>(consts->D, consts->D);
	cstats.log_det_cov = zeros<vec>(1);
	uword cns = 0;
	double clik = 0.0;
	double cprod = 1.0;

	double setlik = 0.0;
	vec likelihood = zeros<vec>(2);
	vec prob_i = zeros<vec>(2);
	uvec n_S = ones<uvec>(2);
	uvec m_S = zeros<uvec>(2);

	int from_class, to_class;

	ivec source_ids_temp = ivec(source_id->begin(), source_id->size());

	// Pick two individuals at random
	uvec inds = zeros<uvec>(2);
	while (inds(0) == inds(1)) {
		inds = randi<uvec>(2, distr_param(0, consts->N - 1));
	}
	// TODO: REMOVE!!!
//	inds << 114 << 67 << endr; // Test case for merging
//	inds << 2 << 3 << endr; // Test case for splitting

	uword s1 = source_ids_temp(inds(0));
	uword s2 = source_ids_temp(inds(1));

	bool merging = (s1 != s2);

	uword n1 = ((uvec) find(source_ids_temp == s1)).size(); 
	uword n2 = ((uvec) find(source_ids_temp == s2)).size();

	// Always merge the smaller class into the largest class
	if (n2 > n1) {
		inds.swap_rows(0, 1);
		uword st = s2;
		s2 = s1;
		s1 = st;
		uword nt = n2;
		n2 = n1;
		n1 = nt;
	}
	D("[split_merge] inds: " << trans(inds));
	D("[split_merge] classes: " << s1 << " " << s2 << "\n");
	D("[split_merge] n: " << n1 << " " << n2 << "\n");
 
	bool s1_baseline = (s1 < consts->nsources);
	bool s2_baseline = (s2 < consts->nsources);
	bool s1_larger = n1 >= n2;
	D("[split_merge] s1 (class " << s1 << ") is " << (s1_baseline ? "" : "NOT ") << "on baseline" << "\n");
	D("[split_merge] s2 (class " << s2 << ") is " << (s2_baseline ? "" : "NOT ") << "on baseline" << "\n");

	D("[split_merge] s" << (s1_larger ? "1" : "2") << " (class " << (s1_larger ? s1 : s2) << ") is larger\n");

	if (merging) {
		/*
		 * 
		 * We're merging. Individuals in source class will be moved to the target
		 * class. If the source class isn't on the baseline it will be deleted and
		 * everyone in the classes above it will be moved back one class
		 *
		 */
		D("[split_merge] Merging\n");

		if (s1_larger) {
			if (!s1_baseline & s2_baseline) {
				from_class = s1;
				to_class = s2;
			} else {
				from_class = s2;
				to_class = s1;
			}
		} else {
			if (!s2_baseline & s1_baseline) {
				from_class = s2;
				to_class = s1;
			} else {
				from_class = s1;
				to_class = s2;
			}
		}
		D("[split_merge]  will consider moving inds from class " << from_class << " to class " << to_class << "\n");
		// Move from_class --> to_class
		source_ids_temp(find(source_ids_temp == from_class)).fill(to_class);
		// If the removed class (from_class) was in the baseline, move the rest of
		// the classes back by one ("compact")
		if (from_class > consts->nsources) {
			source_ids_temp(find(source_ids_temp >= from_class)) -= 1;
		}

	} else {
		/*
		 *
		 * We're splitting. Individual inds[0] will be moved, seeding a new class
		 * (K_plus + 1)
		 *
		 */
		D("[split_merge] Splitting\n");
		from_class = s1;
		to_class = *K_plus + 1;

		// Move the seed individual to the new class
		source_ids_temp[inds[1]] = to_class;

		// Needed for code further down the line
		s2_baseline = false;
	}

	// Get individual and combined stats from the baseline sources
	int a = 0;
	cns = 0;
	if (s1_baseline) {
		n_S(a) += consts->ns(s1);
		sstats[a].means = consts->orgmeans.col(s1).t();
		sstats[a].sum_squares = consts->orgsum_squares.slice(s1);
		
		cns += consts->ns(s1);
		cstats.means += consts->orgmeans.col(s1).t() * consts->ns(s1);
		cstats.sum_squares += consts->orgsum_squares.slice(s1).t();
		++a;
	}
	if (s2_baseline) {
		n_S(a) += consts->ns(s2);
		sstats[a].means = consts->orgmeans.col(s2).t();
		sstats[a].sum_squares = consts->orgsum_squares.slice(s2);

		cns += consts->ns(s2);
		cstats.means += consts->orgmeans.col(s2).t() * consts->ns(s2);
		cstats.sum_squares += consts->orgsum_squares.slice(s2).t();
	}
	if (cns > 0) { cstats.means /= cns; }

	D("[split_merge]   cstats.means: " << cstats.means);
	D("[split_merge]   cstats.sum_squares:\n" << cstats.sum_squares);

	
	// Find all individuals in the classes of inds[0] and inds[1], and shuffle all
	// but the first two
	uvec kk;
	if (merging) {
		kk = shuffle(find((*source_id) == s1 || (*source_id) == s2));
// TODO: REMOVE
//		kk = find((*source_id) == s1 || (*source_id) == s2);
	} else {
		kk = shuffle(find((*source_id) == s2));
// TODO: REMOVE
//		kk = find((*source_id) == s2);
	}
	// Get the first of each class back to the start of the vector
	kk.shed_row(as_scalar(find(kk == inds[0], 1)));
	kk.shed_row(as_scalar(find(kk == inds[1], 1)));
	kk.insert_rows(0, inds);
		
	// 1. pick one at random
	for (uvec::iterator k	= kk.begin(); k != kk.end(); ++k) {
		rowvec y_k = consts->datas.row(*k);

		uword src = ((*source_id)[*k] == s2);

		// 2. update the combined likelihood, counts and stats
		clik += getlik(consts, priors, &cstats, y_k, cns, true, true);
		cns += 1;
//		D("[split_merge] cns: " << cns << "  clik: " << clik);
		update_stats(&cstats, y_k, cns, 1);

		// 3. calculate individual set likelihoods (only for the first two) and skip
		// the rest
		if (*k == kk(0) || *k == kk(1)) {
			int o = (*k == kk(1));
			setlik += getlik(consts, priors, &sstats[o], y_k, n_S[o] - 1, true, true);
			update_stats(&sstats[o], y_k, n_S[o], 1);
			continue;
//			D("\n[split_merge] o: " << o << "  setlik: " << setlik);
		}

//		D("[split_merge] n_S:" << trans(n_S));
		likelihood(0) = getlik(consts, priors, &sstats[0], y_k, n_S[0], true, true);
		likelihood(1) = getlik(consts, priors, &sstats[1], y_k, n_S[1], true, true);

//		D("[split_merge]    lik: " << trans(likelihood));
		vec likelihoods = exp(likelihood);
//		D("[split_merge]   liks: " << trans(likelihoods));

		m_S[0] = s1_baseline ? n_S[0] - consts->ns[s1] : n_S[0];
		m_S[1] = s2_baseline ? n_S[1] - consts->ns[s2] : n_S[1];
		prob_i[0] = (m_S[0] * likelihoods[0]) / sum(m_S % likelihoods);
		prob_i[1] = 1 - prob_i[0];
		
		// If merging, s = src, else take a random pick between 0 and 1
		int s = merging ? src : (as_scalar(randu(1)) > prob_i[0]);

		setlik += likelihood[s];
		cprod *= prob_i[s];
		n_S[s]++;
		update_stats(&sstats[s], y_k, n_S[s], 1);

		if ((!merging) && s == 1) {
			source_ids_temp[*k] = to_class;
		}
	}

	double MH_prior; 
	double MH_lik;
	double MH_rat;

	m_S[0] = s1_baseline ? n_S[0] - consts->ns[s1] : n_S[0];
	m_S[1] = s2_baseline ? n_S[1] - consts->ns[s2] : n_S[1];

	if (merging) {
		MH_prior = exp(lgamma(n1 + n2) - (lgamma(m_S[0]) + lgamma(m_S[1]))) / alpha;
		MH_lik = exp(clik - setlik);
		MH_rat = MH_prior * MH_lik * cprod;
	} else {
		MH_prior = exp(lgamma(m_S[0]) + lgamma(m_S[1]) - lgamma(n2)) * alpha;
		MH_lik = exp(setlik - clik);
		MH_rat = MH_prior * MH_lik / cprod;
	}
	D("[split_merge] * alpha: " << alpha << "  MH_prior: " << MH_prior << "  MH_lik: " << MH_lik << "  MH_rat: " << MH_rat << "\n");

	if (as_scalar(randu(1)) < MH_rat) {
// TODO: REMOVE!!!
//	if (0.5 < MH_rat) {
		D("[split_merge]  really " << (merging ? "merging" : "splitting") << "\n");
		if (merging) {
			// First update suff stats of new merged group
			(*counts)[s2] = n1 + n2;
			(*allcounts)[s2] += n_S[1];
			getlik(consts, priors, &cstats, zeros<rowvec>(consts->D), n1 + n2, false, true);
			(*stats)[s2] = cstats;

			// Then delete the old table
			*K_plus -= 1;

			(*counts).shed_row(s1);
			(*allcounts).shed_row(s1);
			//
			// DROP item s1 from each item on stats and then
			// INSERT new zeroed ones at K_plus + 1
			NORMslice new_NORM;
			new_NORM.means = zeros<rowvec>(consts->D);
			new_NORM.sum_squares = zeros<mat>(consts->D, consts->D);
			new_NORM.inv_cov = zeros<mat>(consts->D, consts->D);
			new_NORM.log_det_cov = zeros<vec>(1);

			(*stats).erase((*stats).begin() + s1);
			(*stats).push_back(new_NORM);

//			D("[split_merge] shedded row: " << trans((*counts).rows(0, *K_plus + 1)) << "\n");
			(*counts).insert_rows(*K_plus + 1, 1, true);
			(*allcounts).insert_rows(*K_plus + 1, 1, true);
//			D("[split_merge] inserted row at " << *K_plus + 1 << ": " << trans((*counts).rows(0, *K_plus + 1)) << "\n");
		} else {
			// Splitting
			D("[split_merge] ** splitting class " << from_class << " into classes " << from_class << " and " << to_class << "\n");
			(*counts)[s2] = n_S[0] - consts->ns[s1] * s1_baseline;
			(*counts)[to_class] = n_S[1];
			(*allcounts)[s2] = n_S[0];
			(*allcounts)[to_class] = n_S[1];
			
			getlik(consts, priors, &sstats[1], zeros<rowvec>(consts->D), (*counts)[to_class], false, true);
			(*stats)[to_class] = sstats[1];
			
			getlik(consts, priors, &sstats[0], zeros<rowvec>(consts->D), (*counts)[s2], false, true);
			(*stats)[s2] = sstats[0];

			*K_plus += 1;
		}	

		*source_id = source_ids_temp;
	} else {
		D("[split_merge]   won't really " << (merging ? "merge" : "split") << "\n");
	}

	D("[split_merge] n[0]: " << (*counts)[0] << "  means: " << (*stats)[0].means);
	D("[split_merge] n[1]: " << (*counts)[1] << "  means: " << (*stats)[1].means);
	D("[split_merge] n[2]: " << (*counts)[2] << "  means: " << (*stats)[2].means);
	D("[split_merge] n[3]: " << (*counts)[3] << "  means: " << (*stats)[3].means);
}

