/*
 * onepl.cpp
 *
 *  Created on: Sep 20, 2016
 *      Author: milder
 */

#include "onepl.h"

namespace latentregpp {

namespace dichotomous {

onepl::onepl() : model(model_type::onepl, ONE_PARAMETER) {}

double onepl::P(std::vector<double> &theta, const optimizer_vector &parameters, int i ) {
	//Initialized with gamma_k value
	double eta = parameters(0);

	//Computing dot product
	for ( size_t j = 0; j < theta.size(); ++j )
		eta += 1 * theta[j]; //no alpha in this model

	double P = 1.0 / (1.0 + std::exp(-eta));
	P = std::max(P, LOWER_BOUND_);
	P = std::min(P, UPPER_BOUND_);
	return P;
}

onepl::~onepl() {}

}

} /* namespace latentregpp */
