/*
 * mstep.cpp
 *
 *  Created on: 13/04/2016
 *      Author: Milder
 */

#include "../../dichotomous/estimation/mstep.h"

namespace latentregpp {

namespace dichotomous {

Qi::Qi (int i, estimation_data *d) : i(i), data(d) { }

double Qi::operator() ( const optimizer_vector& item_i ) const {
	//Value of Qi
	double value = 0;
	//Number of quadrature points
	int G = data->G;
	//Matrix r
	matrix<double> &r = data->r;
	//Model used
	model &m = data->m;
	//Latent trait vectors
	matrix<double> &theta = data->theta;
	//f
	std::vector<double> &f = data->f;

	for ( int g = 0; g < G; ++g ) {
		std::vector<double> &theta_g = *theta.get_pointer_row(g);
		double P_gi = m.P(theta_g, item_i);
		value += r(g, i) * log(P_gi) + (f[g] - r(g, i)) * log(1 - P_gi);
	}

	return value;
}

double Mstep(estimation_data &data, int current) {
	double max_difference = 0.0;
	int next = (current + 1) % ACCELERATION_PERIOD;

	int &p = data.p;
	model &m = data.m;
	std::vector<optimizer_vector> &current_zeta = data.zeta[current];
	std::vector<optimizer_vector> &next_zeta = data.zeta[next];

	if ( next_zeta.empty() ) {
		for ( auto c : current_zeta )
			next_zeta.push_back(c);
	}

	std::set<int> &pinned_items = data.pinned_items;

	#pragma omp parallel for schedule(dynamic) reduction(max:max_difference)
	for ( int i = 0; i < p; ++i ) {
		/**
		 * If it is multidimensional and this is one of the pinned items
		 * i.e the first item of a dimension
		 * this item is just skipped
		 * */

		if ( pinned_items.count(i) ) continue;


		next_zeta[i] = current_zeta[i];

		//Calling BFGS from dlib to optimize Qi with approximate derivatives
		dlib::find_max_using_approximate_derivatives(dlib::bfgs_search_strategy(),
					   dlib::objective_delta_stop_strategy(OPTIMIZER_DELTA_STOP),
					   Qi(i, &data), next_zeta[i], -1);

		//Computing difference of current item
		if ( m.parameters < 3 ) {
			for ( int j = 0; j < next_zeta[i].size(); ++j )
				max_difference = std::max(max_difference, std::abs(next_zeta[i](j) - current_zeta[i](j)));
		} else {
			for ( int j = 0; j < next_zeta[i].size() - 1; ++j )
				max_difference = std::max(max_difference, std::abs(next_zeta[i](j) - current_zeta[i](j)));
			double c_current = current_zeta[i](current_zeta[i].size() - 1);
			double c_next = next_zeta[i](next_zeta[i].size() - 1);

			c_current = 1.0 / (1.0 + exp(-c_current));
			c_next = 1.0 / (1.0 + exp(-c_next));

			max_difference = std::max(max_difference, std::abs(c_next - c_current));
		}
	}

	return max_difference;
}

}

} /* namespace latentregpp */

