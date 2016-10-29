/*
 * threepl.cpp
 *
 *  Created on: Sep 20, 2016
 *      Author: milder
 */

#include "threepl.h"

namespace latentregpp {

namespace dichotomous {

threepl::threepl() : model(model_type::threepl, THREE_PARAMETERS) {}

double threepl::P(std::vector<double> &theta, const optimizer_vector &parameters, int i) {
	double gamma_parameter = parameters(parameters.size() - 1);

	//uncommented line below for reparameter a c value [0,1] from gamma in R
	double c = 1.0 / (1.0 + exp(-gamma_parameter));

	//Initialized with gamma value
	double eta = parameters(parameters.size() - 2);

	//Computing dot product
	for ( size_t j = 0; j < theta.size(); ++j )
		eta += parameters(j) * theta[j];

	/**three different formulas**/

	//with clear
	//double P = (1.0 / (1.0 + exp(-gamma_parameter))) + (1.0 - (1.0 / (1.0 + exp(-gamma_parameter)))) / (1.0 + exp(-eta));

	//without clear
	//double P = (1.0 / (1.0 + std::exp(-gamma_parameter))) + (1.0 / (1.0 + std::exp(gamma_parameter))) * (1.0 / (1.0 + std::exp(-eta)));

	//with reparameter
	double P = c + (1.0 - c) / (1.0 + exp(-eta));

	P = std::max(P, LOWER_BOUND_);
	P = std::min(P, UPPER_BOUND_);
	return P;
}

threepl::~threepl() {}

}

} /* namespace latentregpp */
