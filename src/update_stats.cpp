#include "unmixR.h"

using namespace arma;

// CHECKED

void update_stats(NORMslice *stats, rowvec y, int counts, int m) {
	if (m < 0) {
		stats->means = (1.0 / counts) * ((counts + 1.0) * stats->means - y);
		stats->sum_squares -= trans(y) * y;
	} else {
		stats->means += (1.0 / counts) * (y - stats->means);
		stats->sum_squares += trans(y) * y;
	}
}

