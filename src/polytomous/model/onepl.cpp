/*
 * onepl.cpp
 *
 *  Created on: Sep 20, 2016
 *      Author: milder
 */

#include "onepl.h"

namespace latentregpp {

namespace polytomous {

onepl::onepl(int d, std::vector<int> *categories_item) :
	model(model_type::onepl, ONE_PARAMETER, d, categories_item) {
}

onepl::~onepl() {}


double onepl::Pstar_ik(std::vector<double> &theta, const optimizer_vector &parameters, int i, int k) {
	int mi = (*categories_item)[i];
	/**
	 * Base cases
	 *
	 * Because of implementation indexes start from 0 and theory from 1
	 * It's necessary subtract 1
	 * */

	if ( k == -1 ) return 1;
	if ( k == mi - 1 ) return 0;

	//Initialized with gamma_k value
	double eta = parameters(k);

	//Computing dot product
	for ( size_t i = 0; i < theta.size(); ++i )
		eta += 1 * theta[i]; //no alpha in this model

	return 1.0 / (1.0 + std::exp(-eta));
}

double onepl::Pik(std::vector<double> &theta, const optimizer_vector &parameters, int i, int k) {
	double P_ik = Pstar_ik(theta, parameters, i, k - 1) - Pstar_ik(theta, parameters, i, k);

	P_ik = std::max(P_ik, LOWER_BOUND_);
	P_ik = std::min(P_ik, UPPER_BOUND_);
	return P_ik;
}

}

} /* namespace latentregpp */
