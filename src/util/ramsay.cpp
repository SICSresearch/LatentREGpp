/*
 * squarem.cpp
 *
 *  Created on: Jul 28, 2016
 *      Author: Milder
 */

#include "ramsay.h"

namespace latentregpp {

void ramsay(std::vector<optimizer_vector> zeta[ACCELERATION_PERIOD],
		std::set<int> &pinned) {

	double num = 0, den = 0, accel = 0;
	int p = zeta[0].size();

	for ( int i = 0; i < p; ++i ) {
		if ( pinned.count(i) ) continue;
		for ( int j = 0; j < zeta[0][i].size(); ++j ) {
			double dx = zeta[2][i](j) - zeta[1][i](j);
			double dx2 = zeta[1][i](j) - zeta[0][i](j);
			double d2x2 = dx - dx2;

			num += dx * dx;
			den += d2x2 * d2x2;
		}
	}

	accel = std::max(1 - sqrt(num / den), MINIMUM_ACCEL);

	#pragma omp parallel for schedule(dynamic)
	for ( int i = 0; i < p; ++i ) {
		if ( pinned.count(i) ) continue;
		for ( int j = 0; j < zeta[0][i].size(); ++j )
			zeta[2][i](j) = (1 - accel) * zeta[2][i](j) + accel * zeta[1][i](j);
	}
}

} /* namespace latentregpp */
