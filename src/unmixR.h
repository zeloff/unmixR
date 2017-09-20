#include <RcppArmadillo.h>

#define NDEBUG

#ifdef NDEBUG
# define D(X) Rcpp::Rcout << X
#else
# define D(X) do {} while (0)
#endif

struct STUD { // consts
	arma::mat datas;
	arma::mat baseline;
	arma::mat orgmeans;
	arma::cube orgsum_squares;
	arma::cube orginv_cov;
	arma::vec orglog_det_cov;
	arma::uvec ns;
	arma::ivec label;
	Rcpp::CharacterVector labels;
	arma::ivec sources;
	int nsources;
	int N;

	arma::vec pc_max_ind;
	arma::vec pc_gammaln_by2;
	double pc_log_pi;
	arma::vec pc_log;
	int D;
} ;

struct MNIW { // priors
	double a_0;
	double b_0;
	double k_0;
	double ak_0;
	double bk_0;
	double v_0;
	arma::rowvec mu_0;
	arma::mat lambda_0;
} ;

struct NORM { // stats
	arma::mat means;       
	arma::cube sum_squares;
	arma::cube inv_cov;
	arma::vec log_det_cov;
} ;

struct NORMslice {
	arma::rowvec means;
	arma::mat sum_squares;
	arma::mat inv_cov;
	arma::vec log_det_cov;
} ;

void crp_gibbs(arma::ivec*, STUD*, MNIW*, std::vector<NORMslice>*, arma::ivec*, arma::ivec*, int*, double, arma::vec);

double getlik(STUD *consts, MNIW *priors, NORMslice *stats, arma::rowvec y, int n, bool lik, bool suffs);

arma::vec p_for_1(STUD*, MNIW*, int, arma::mat);

arma::mat rwish(double, arma::mat);
SEXP rwishart(double, Rcpp::NumericMatrix);
arma::mat rwishart(double, arma::mat);
SEXP rinvwishart(double, Rcpp::NumericMatrix);
arma::mat rinvwishart(double, arma::mat);
void print_stats(std::vector<NORMslice>*, char[]);

void split_merge(arma::mat, arma::ivec*, STUD*, int, MNIW*, std::vector<NORMslice>*, arma::ivec*, arma::ivec*, int*, double);

void update_stats(NORMslice*, arma::rowvec, int, int);

void update_alpha(double*, int, int, MNIW);

void update_prior(STUD, int, std::vector<NORMslice>, arma::ivec, MNIW*);

void preallocate(STUD*, std::vector<NORMslice>*, MNIW*, arma::ivec *allcounts);
void preclass(STUD*, std::vector<NORMslice>*, MNIW*, arma::ivec*, arma::ivec*, arma::ivec*);


SEXP struct_list(STUD*, std::vector<NORMslice>*, arma::ivec*);

