#include <math.h>
#include <stdlib.h>
#include <iostream>

#include "unmixR.h"

// Comment out this line for debugging output:
# define D(X) do {} while (0)

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
SEXP DPM_sampler(int thin, int iters, NumericMatrix data_, NumericMatrix baseline_, CharacterVector labels_, List consts_, List params) {

	int D;
	int N;

	D("[DPM_sampler] Starting...\n");
	mat data(data_.begin(), data_.nrow(), data_.ncol(), false);

	try {
		D = data.n_cols;
		N = data.n_rows;
	} catch ( std::exception& __ex__ ) {
		std::cerr << "Data matrix has invalid dimensions. Please check inputs.\n";
		forward_exception_to_r( __ex__ );
	}

	STUD consts;
	MNIW priors;
	std::vector<NORMslice> stats(N);

	try {
		consts.datas = data;
		consts.baseline = mat(baseline_.begin(), baseline_.nrow(), baseline_.ncol(), false);

		consts.labels = sort_unique(labels_);
		IntegerVector label_idx = match(labels_, consts.labels) - 1;
		consts.label = ivec(label_idx.begin(), label_idx.size(), false);

		consts.sources = ivec(label_idx.begin(), label_idx.size(), false);
		std::unique(consts.sources.begin(), consts.sources.end());
		consts.sources.resize(consts.labels.size());

		D("[DPM_Sampler] labels: " << consts.labels << "\n");
		D("[DPM_Sampler] label: " << trans(consts.label) << "\n");
		D("[DPM_Sampler] sources: " << trans(consts.sources) << "\n");

		consts.nsources = consts.sources.size();
		consts.N = N;

		/*
		 * set by preallocate()
		 */
		consts.ns.zeros(consts.nsources);
		consts.orgsum_squares.zeros(D, D, consts.nsources);
		consts.orgmeans.zeros(D, consts.nsources);
		consts.orginv_cov.zeros(D, D, consts.nsources);
		consts.orglog_det_cov.zeros(consts.nsources);

		NumericVector pcmaxind = consts_("pc_max_ind");
		consts.pc_max_ind = vec(pcmaxind.begin(), pcmaxind.size(), false);

		NumericVector pcglnb2 = consts_("pc_gammaln_by2");
		consts.pc_gammaln_by2 = vec(pcglnb2.begin(), pcglnb2.size(), false);

		consts.pc_log_pi = log(M_PI); 

		NumericVector pclog = consts_("pc_log");
		consts.pc_log = vec(pclog.begin(), pclog.size(), false);

		consts.D = D;
	} catch ( std::exception& __ex__ ) {
		std::cerr << "Error initializing constants. Please check 'consts'.\n";
		forward_exception_to_r( __ex__ );
	}

	try {
		priors.a_0 = params("a_0");
		priors.b_0 = params("b_0");
		priors.k_0 = params("k_0");
		priors.ak_0 = params("ak_0");
		priors.bk_0 = params("bk_0");
		priors.v_0 = params("v_0");
		NumericVector mu_0_ = params("mu_0");
		priors.mu_0 = rowvec(mu_0_.begin(), mu_0_.size(), false);
		NumericMatrix lambda_0_ = (SEXP) params("lambda_0");
		priors.lambda_0 = mat(lambda_0_.begin(), lambda_0_.nrow(), lambda_0_.ncol(), false);
	} catch ( std::exception& __ex__ ) {
		std::cerr << "Error initializing priors. Please check 'params'.\n";
		forward_exception_to_r( __ex__ );
	}

	for (int n = 0; n < N; n++) {
		stats[n].means = rowvec(D);
		stats[n].means.zeros();
		stats[n].sum_squares = mat(D, D);
		stats[n].sum_squares.zeros();
		stats[n].inv_cov = mat(D, D);
		stats[n].inv_cov.zeros();
		stats[n].log_det_cov = vec(1);
		stats[n].log_det_cov.zeros();
	}


	double alpha = 0.1;

	// "initialize structures"

	uword i;
	uword nit = floor(iters / thin + 0.5); // both are positive ints, so this is ok for rounding to int
  bool drop;
	
	uword trace;
	uvec iter_trace;
	iter_trace.zeros(nit);

	vec k_0s(nit);
	imat class_ids(consts.N, nit);
	ivec K_record(nit);
	vec alpha_record(nit);

	iters--;
	int max_class_id = consts.N;
	ivec class_id(consts.N);
	class_id.fill(-1);
	class_id[0] = 0;
	ivec counts = zeros<ivec>(max_class_id);     // Counts excluding baseline
	ivec allcounts = zeros<ivec>(max_class_id);  // Counts including baseline
	int K_plus = consts.nsources - 1;
	vec p_under_prior_alone(consts.N);

	preallocate(&consts, &stats, &priors, &allcounts);
	D("[main] preallocate: " << trans(counts.head_rows(5)));
	D("[main]   * K_plus " << K_plus << "\n");
	
	for (int i=0; i < consts.nsources; ++i) {
		D("[main] " << i << "  source: " << consts.sources(i));
		D("  label: " << consts.labels[i] << "\t");
		D("  ns: " << consts.ns[i]);
		D("  orgmeans: " << trans(consts.orgmeans.col(i)));
	}

	preclass(&consts, &stats, &priors, &allcounts, &counts, &class_id);

	D("[main] preclass: " << trans(counts.head_rows(5)));
	D("[main]   * K_plus " << K_plus << "\n");

	for (int iter = 0; iter <= iters; iter++) {
		D("\n[main] starting iteration " << iter << "\n");
		D("[main] calculating p under prior alone\n");

    drop = false;
    trace = 0;
	  
		try {
			p_under_prior_alone = p_for_1(&consts, &priors, N, data);

			D("[main] calling the Gibbs sampler\n");
			D("[main]   * K_plus " << K_plus << "  counts: " << trans(counts.head(4)));
      trace = 1;
			crp_gibbs(&class_id, &consts, &priors, &stats, &counts, &allcounts, &K_plus, alpha, p_under_prior_alone);

			D("[main] Gibbs sampler is done.\n");
			D("[main] gibbs: " << trans(counts.head_rows(5)));

			D("[main] Splitting / Merging.\n");
			D("[main]   * K_plus " << K_plus << "  counts: " << trans(counts.head(4)));
      trace = 2;
			try {
				split_merge(data, &class_id, &consts, N, &priors, &stats, &counts, &allcounts, &K_plus, alpha);
			} catch ( std::exception& __ex__ ) {
				std::cerr << "Error splitting/merging.\n";
				forward_exception_to_r( __ex__ );
			}

			D("[main] Splitting / Merging is done.\n");

			D("[main]   *         K_plus: " << K_plus << "\n");
			D("[main]   *         counts: " << trans(counts.head(4)));
			D("[main]   *      allcounts: " << trans(allcounts.head(4)));
			D("[main]   * stats.means[0]: " << stats[0].means);
			D("[main]   * stats.means[1]: " << stats[1].means);
			D("[main]   * stats.means[2]: " << stats[2].means);
			D("[main]   * stats.means[3]: " << stats[3].means);
			D("[main]   *    priors: k_0: " << priors.k_0 << "  mu_0: " << priors.mu_0);

			D("[main] Updating alpha.\n");
      trace = 3;
			update_alpha(&alpha, N, K_plus, priors);

			D("[main] Alpha updated.\n");
			D("[main]   * alpha: " << alpha << "\n");

			trace = 4;

			D("[main] Updating prior.\n");
			update_prior(consts, K_plus, stats, counts, &priors);
			D("[main]   * k_0: " << priors.k_0 << "  mu_0: " << priors.mu_0);
		} catch (std::domain_error& de) {
		  std::cerr << "[main] Domain error at step " << trace << " of iteration " << i << ".\n";
		  forward_exception_to_r(de);
		} catch ( std::exception& __ex__ ) {
			std::cerr << "[main] Error at step " << trace << " of iteration " << i << ".\n";
		  drop = true;
		}

		i = floor(iter / thin + 0.5);
    
//		if (iter % 1000 == 0) {
			D("[main] Iteration " << iter << " (" << i << "/" << nit << ") done");
//		};

		if (!drop) {
		  iter_trace(i) = i;  
		}
		k_0s(i) = priors.k_0;
		class_ids.col(i) = class_id;
		K_record(i) = K_plus + 1;
		alpha_record(i) = alpha;
		D(".\n");
	}

	D("[main] All done.\n");
	
	iter_trace(0) = 1;
	uvec iter_trace_ = find(iter_trace > 0);
	iter_trace_(0) = 0;
	
	vec k0s(k_0s.elem(iter_trace_));
	imat ids(class_ids.cols(iter_trace_));
	ivec ks(K_record.elem(iter_trace_));
	vec alphas(alpha_record.elem(iter_trace_));
	
	D("[main] Dropped " << nit - iter_trace_.n_rows << " iterations.\n");
	/*
	List lstats = List::create(
			Named("means") = xs_means,
			Named("sum_squares") = xs_ss,
			Named("inv_cov") = xs_inv_cov,
			Named("log_det_cov") = xs_ldc);
			*/
	List ret = List::create(
	    Named("iters") = iter_trace_ + 1,
			Named("class_ids") = ids,
			Named("K_record") = ks,
			Named("k_0s") = k0s + 1,
			Named("alphas") = alphas);

	return(ret);
}

