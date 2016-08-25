/*
 * pi_matrix_test.cpp
 *
 *  Created on: 25/04/2016
 *      Author: Milder
 */

#include "test.h"
#include <iostream>

namespace latentregpp {

	bool test_pi ( matrix<double> &pi ) {
		const double eps = 1e-5;

		bool correct = true;
		for ( int j = 0; j < pi.columns(0); ++j ) {
			double sum = 0;
			for ( int i = 0; i < pi.rows(); ++i )
				sum += pi(i, j);
			correct &= 100 == int((sum + eps) * 100);
		}

		return correct;
	}

	bool test_r ( std::vector<matrix<double> > &r, int N, int p ) {
		double sum = 0;
		for ( size_t i = 0; i < r.size(); ++i )
			for ( int j = 0; j < r[i].rows(); ++j )
				for ( int k = 0; k < r[i].columns(j); ++k )
					sum += r[i](j, k);
		const double eps = 1e-5;
		return N * p == int(sum + eps);
	}

	bool test_r ( matrix<double> &r, int N, int p ) {
		double sum = 0;
		for ( int j = 0; j < r.rows(); ++j )
			for ( int k = 0; k < r.columns(j); ++k )
					sum += r(j, k);
		const double eps = 1e-5;

		return N * p == int(sum + eps);
	}

}



