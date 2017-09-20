#include "unmixR.h"

using namespace Rcpp;
using namespace arma;

void crp_gibbs(ivec* source_id, STUD *consts, MNIW *priors, std::vector<NORMslice>* stats, ivec *counts, ivec *allcounts, int *K_plus, double alpha, vec pupa)
{
	int trace = 0;

	try {

		for (int i = 0; i < consts->N; i++) {

			trace = 1;

			rowvec y = consts->datas.row(i);

			int old_source_id = as_scalar(source_id->at(i));

			trace = 2;

			// Temporary variables
			std::vector<NORMslice> tempstats(consts->N);
			for (int n = 0; n < consts->N; n++) {
				tempstats[n].means = (*stats)[n].means;
				tempstats[n].sum_squares = (*stats)[n].sum_squares;
				tempstats[n].inv_cov = (*stats)[n].inv_cov;
				tempstats[n].log_det_cov = (*stats)[n].log_det_cov;
			}

			trace = 3;

			int last_class = *K_plus;
			ivec source_id_temp(*source_id);
			ivec tempcounts(*counts);
			ivec tempallcounts(*allcounts);

			bool dels = false;

			trace = 4;

			// Remove this individual from that source's count
			tempcounts(old_source_id)--; 
			tempallcounts(old_source_id)--; 

			if ((tempcounts[old_source_id] == 0) & (old_source_id > (consts->nsources - 1))) { // If this source is now empty and is the last source
				dels = true; // then something will change
				//
				// Drop this source and move everyone back (no zeros are allowed in the
				// middle of the tempcounts and tempallcounts vectors)
				tempcounts.shed_row(old_source_id);
				tempcounts.insert_rows(last_class, 1, true);
				tempallcounts.shed_row(old_source_id);
				tempallcounts.insert_rows(last_class, 1, true);

				NORMslice new_NORM;
				new_NORM.means = zeros<rowvec>(consts->D);
				new_NORM.sum_squares = zeros<mat>(consts->D, consts->D);
				new_NORM.inv_cov = zeros<mat>(consts->D, consts->D);
				new_NORM.log_det_cov = zeros<vec>(1);

				tempstats.erase(tempstats.begin() + old_source_id);
				tempstats.push_back(new_NORM);

				// Find out how many individuals were in this source and in the ones
				// above it and move them to the source below them
				source_id_temp.elem(find(source_id_temp >= old_source_id)) -= 1;

				// "Shorten" the number of total sources, since we moved back all
				// individuals on the top sources
				last_class--;
				//
				// The old source is now the top source
				old_source_id = last_class + 1;
			} else { // This source wasn't empty
				// Update stats
				update_stats(&tempstats[old_source_id], y, tempallcounts(old_source_id), -1);
			} // if

			trace = 5;

			/*
			 * Prior with new source probability
			 */
			/* Complete the new CRP prior with the new source prob */
			vec tt(conv_to<vec>::from(tempcounts.rows(0, last_class + 1)));
			uvec tt_zeros = find(tt == 0);
			double a_ = alpha / (tt_zeros.n_elem + 1);
			tt.elem(tt_zeros).fill(a_);
			tt(tt.n_elem - 1) = a_;

			vec prior = tt / (consts->N - 1 + alpha);

			trace = 6;

			/*
			 * New likelihood
			 */
			vec likelihood = zeros(prior.size());
			for (int ell = 0; ell <= last_class; ell++) {
//				D("[crp_gibbs] i: " << i << "  ell: " << ell << "  old_source_id: " << old_source_id << " - " << (old_source_id == ell || old_source_id == -1) << "\n");
				likelihood[ell] = getlik(consts, priors, &tempstats[ell], y, tempallcounts[ell], true, (old_source_id == ell || old_source_id == -1));
			}
			likelihood(last_class + 1) = pupa(i); // actually, this is new 'top source' which was just added
			likelihood = exp(likelihood - max(likelihood));
			likelihood /= sum(likelihood);

			trace = 7;

			/*
			 * Posterior
			 */
			vec posterior(prior.size());
			posterior = prior % likelihood;
			posterior /= sum(posterior); // Normalize the posterior

			D("[crp_gibbs] posterior: " << trans(posterior));
			/*
			 * Pick the new source
			 */

			trace = 8;

			vec cs = cumsum(posterior);
			double r = as_scalar(randu(1));
			// TODO: remove!!!
//			double r = 0.5;
			D("[crp_gibbs] cs: " << trans(cs) << "  r: " << r << "\n");
			uvec ff(find(cs > r));
			int new_source_ids = ff(0);

			trace = 9;

			tempcounts(new_source_ids)++;
			tempallcounts(new_source_ids)++;
			bool newc = false;
			if (new_source_ids == last_class + 1) {
				newc = true;
				last_class++;
			}

			/*
			 * If the new source_id != old source_id without any rearrangements or if the
			 * new source_id != K_plus+1 with rearrangement, update things. Else, do nothing
			 */
			if (((old_source_id != new_source_ids) & !dels) | (!newc & dels)) {
				trace = 10;
				// Store the new id
				source_id_temp(i) = new_source_ids;

				// Temporary values are now definitive
				*source_id = source_id_temp;
				*K_plus = last_class;
				*counts = tempcounts;
				*allcounts = tempallcounts;

				trace = 11;
				update_stats(&tempstats[new_source_ids], y, (*allcounts)(new_source_ids), 1);
				trace = 12;
				getlik(consts, priors, &tempstats[new_source_ids], y, (*allcounts)(new_source_ids), false, true);
				trace = 13;
				for (int ts = 0; ts < consts->D; ts++) {
					(*stats)[ts].means = tempstats[ts].means;
					(*stats)[ts].sum_squares = tempstats[ts].sum_squares;
					(*stats)[ts].inv_cov = tempstats[ts].inv_cov;
					(*stats)[ts].log_det_cov = tempstats[ts].log_det_cov;
				}
			} // if changed 
			trace = 14;
		} // for
		D("\n");
	} catch ( std::domain_error& __de__ ) {
		D("[crp_gibbs] error after step " << trace << "\n");
		throw __de__;
	} catch ( std::exception& __ex__ ) {
		D("[crp_gibbs] error after step " << trace << "\n");
		throw __ex__;
	}
}
