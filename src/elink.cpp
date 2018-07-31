#include "unmixR.h"

// Comment out this line for debugging output:
# define D(X) do {} while (0)

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
SEXP elink(NumericMatrix ss)
{
	mat S(ss.begin(), ss.nrow(), ss.ncol(), false);

	uword n = S.n_rows;
	uword niter = S.n_cols;

	D("[elink] n: " << n << "\n");
	D("[elink] niter: " << niter << "\n");

	mat meanDI;

	umat Z;		// "Best pairs"...
	vec z;    // ... and its distance

	umat R;		// Leafs and groups
	uvec MID;	// Current groups
	
	Z.zeros(n - 1, 2);
	z.zeros(n - 1);

	R.zeros(n + n - 1, n - 1);
	IntegerVector h1 = Range(0, n - 1);
	R.submat(0, 0, n - 1, 0) = as<uvec>(h1);

	meanDI.zeros(n, n);
	MID = as<uvec>(h1);

	double sum_, max_ = 0;
	for (uword i = 0; i < n - 1; i++) {
		for (uword j = i + 1; j < n ; j++) {
			sum_ = sum(S.row(i) == S.row(j)) / (niter * 1.0);
			meanDI(i, j) = sum_;
			meanDI(j, i) = sum_;
			if (sum_ > max_) {
				max_ = sum_;
				Z(0, 0) = i;
				Z(0, 1) = j;
				z(0) = 1 - sum_;
			}
		}
	}

	D("[elink] First pair: " << Z(0, 0) << " " << Z(0, 1) << " " << z(0) << "\n");
	// New group
	R.submat(n, 0, n, 1) = Z.row(0);

	// Drop best pair from distance matrix and from MID
	// j > i, so Z(1) > Z(0), so Z(1) must be taken out first
	meanDI.shed_row(Z(1));
	meanDI.shed_row(Z(0));
	meanDI.shed_col(Z(1));
	meanDI.shed_col(Z(0));

	MID.shed_row(Z(1));
	MID.shed_row(Z(0));

	//
	
	uvec an0, an1, allnode;

	for (uword s = 1; s < n - 1; s++) {
		D("[elink] s: " << s << "\n");
		D("[elink]   Z(s - 1): " << Z.row(s - 1));
		// Get the nodes on the same group as the ones on Z.row(s - 1) 
		an0 = trans(R.row(Z(s - 1, 0)));
		an0 = an0.elem(find(an0));
		an1 = trans(R.row(Z(s - 1, 1)));
		an1 = an1.elem(find(an1));
		allnode = join_cols(an0, an1);
		D("[elink]   Nodes on the same group as Z(s - 1): " << trans(allnode));

		// Add new group n + s - 1
		D("[elink]   New group: " << n + s - 1 << "\n");
		MID.insert_rows(MID.n_rows, 1);
		MID(MID.n_rows - 1) = n + s - 1;
		meanDI.insert_rows(meanDI.n_rows, 1, true);
		meanDI.insert_cols(meanDI.n_cols, 1, true);
		vec DIC = ones<vec>(n - s - 1);

		// Calculate the distance of each leaf/group to the newly added group
		D("[elink]   Calculating distance to other leafs/groups (0 - " << n - s - 2 << ")\n");
		uvec nodes;
		double temp;
		for (uword o = 0; o < n - s - 1; o++) {
			nodes = trans(R.row(MID(o)));
			nodes = nodes.elem(find(nodes));
			for (uword y = 0; y < allnode.n_rows; y++) {
				for (uword z = 0; z < nodes.n_rows; z++) {
					temp = sum(S.row(allnode(y)) == S.row(nodes(z))) / (niter * 1.0);
					if (temp < DIC(o)) {
						DIC(o) = temp;
					}
				}
			}
		}
		D("[elink]   DIC(" << DIC.n_rows << "): " << sum(DIC) << "\n");
		D("[elink]   Rows: 0:" << n - s - 2 << ",\tCol: " << n - s - 1 << "\n");
		D("[elink]    Row: " << n - s - 1 << ",\tCol: 0:" << n - s - 2 << "\n");

		meanDI.submat(0, n - s - 1, n - s - 2, n - s - 1) = DIC;
		meanDI.submat(n - s - 1, 0, n - s - 1, n - s - 2) = trans(DIC);

		// Find new closest pair...
		D("[elink]   Finding closest pair (ranges: 0-" << n - s - 2 << ", i-" << n - s - 1 << "): ");
		int ro = 0, col = 0;
		for (uword i = 0; i < n - s - 1; i++) {
			for (uword j = i; j < n - s; j++) {
				if (meanDI(i, j) >= meanDI(ro, col)) {
					ro = i;
					col = j;
				}
			}
		}
		D(ro << "-" << col << " (" << meanDI(ro, col) << ")\n");

		// ... store it ...
		Z(s, 0) = MID(ro);
		Z(s, 1) = MID(col);
		z(s) = 1.0 - meanDI(ro, col);

		// ... get the elements in the same group ...
		uvec newgr1, newgr2;
		newgr1 = trans(R.row(MID(ro)));
		newgr1 = newgr1.elem(find(newgr1));
		newgr2 = trans(R.row(MID(col)));
		newgr2 = newgr2.elem(find(newgr2));

		if (s == n - 2) {
			break;
		}

		// ... update the R matrix ...
//		D("[elink]   newgr1: " << trans(newgr1) << "[elink]   newgr2: " << trans(newgr2));
		if (newgr1.n_rows + newgr2.n_rows > 1) {
			R.submat(n + s, 0, n + s, newgr1.n_rows + newgr2.n_rows - 1) = trans(join_cols(newgr1, newgr2));
		}
		D("[elink]  R(" << n + s << ", 0:" << newgr1.n_rows + newgr2.n_rows - 1 << "):" << R.submat(n + s, 0, n + s, newgr1.n_rows + newgr2.n_rows - 1));

		D("[elink]  Dropping " << ro << " and " << col << "\n");
		// ... and finally drop the new best pair
		if (ro > col) {
			meanDI.shed_row(ro);
			meanDI.shed_row(col);
			meanDI.shed_col(ro);
			meanDI.shed_col(col);
			MID.shed_row(ro);
			MID.shed_row(col);
		} else {
			meanDI.shed_row(col);
			meanDI.shed_row(ro);
			meanDI.shed_col(col);
			meanDI.shed_col(ro);
			MID.shed_row(col);
			MID.shed_row(ro);
		}
	}

	List ret = List::create(
			Named("pairs") = Z,
			Named("distances") = z);
	return wrap(ret);
}

