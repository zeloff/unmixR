#include "unmixR.h"

using namespace Rcpp;

/*
SEXP struct_list(STUD *consts, std::vector<NORMslice> *stats, arma::ivec *class_ids) {
	List const_list = List::create(
			Named("datas") = (*consts).datas,
			Named("baseline") = (*consts).baseline,
			Named("orgsum_squares") = (*consts).orgsum_squares,
			Named("orgmeans") = (*consts).orgmeans,
			Named("orginv_cov") = (*consts).orginv_cov,
			Named("orglog_det_cov") = (*consts).orglog_det_cov,
			Named("ns") = (*consts).ns,
			Named("labels") = (*consts).labels,
			Named("sources") = (*consts).sources,
			Named("nsources") = (*consts).nsources,
			Named("N") = (*consts).N,
			Named("D") = (*consts).D,
			Named("pc_max_ind") = (*consts).pc_max_ind,
			Named("pc_gammaln_by2") = (*consts).pc_gammaln_by2,
			Named("pc_log_pi") = (*consts).pc_log_pi,
			Named("pc_log") = (*consts).pc_log
			);

	List stats_list = List::create();
	for (std::vector<NORMslice>::iterator it = (*stats).begin(); it != (*stats).end(); ++it) {
		List sl = List::create(
				Named("means") = (*it).means,
				Named("sum_squares") = (*it).sum_squares,
				Named("inv_cov") = (*it).inv_cov,
				Named("log_det_cov") = (*it).log_det_cov
				);
		stats_list.push_back(sl);
	}

	List ret = List::create(
	    Named("consts") = const_list,
			Named("stats") = stats_list,
			Named("class_ids") = (*class_ids)
			);

	return ret;
}
*/
