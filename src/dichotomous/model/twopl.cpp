/*
 * twopl.cpp
 *
 *  Created on: Sep 20, 2016
 *      Author: milder
 */

#include "twopl.h"

namespace latentregpp {

namespace dichotomous {

twopl::twopl() : model(model_type::twopl, TWO_PARAMETERS) {}

double twopl::P(std::vector<double> &theta, const optimizer_vector &parameters, int i) {
	//Initialized with gamma value
	double eta = parameters(parameters.size() - 1);

	//Computing dot product
	for ( size_t j = 0; j < theta.size(); ++j )
		eta += parameters(j) * theta[j];

	double P = 1.0 / (1.0 + std::exp(-eta));
	P = std::max(P, LOWER_BOUND_);
	P = std::min(P, UPPER_BOUND_);
	return P;
}

twopl::~twopl() {}

}

} /* namespace latentregpp */
