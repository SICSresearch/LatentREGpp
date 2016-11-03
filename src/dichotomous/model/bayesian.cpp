/*
 * bayesian.cpp
 *
 *  Created on: 26/10/2016
 *      Author: jhonatan
 */

#include "bayesian.h"

namespace latentregpp {

namespace dichotomous {

//bayesian::bayesian() : model(model_type::bayesian, THREE_PARAMETERS) {}

bayesian::bayesian(matrix<double> &c) : model(model_type::bayesian, THREE_PARAMETERS) {
	this->c_values = c;
  //this->type = model_type::bayesian;
  //this->parameters = THREE_PARAMETERS;
}

double bayesian::P(std::vector<double> &theta, const optimizer_vector &parameters) {
	double gamma_parameter = parameters(parameters.size() - 1);

	//uncommented line below for reparameter a c value [0,1] from gamma in R
	//double c = 1.0 / (1.0 + exp(-gamma_parameter));

	//need to put c initial value here
	int item = 0;
	double c = 0.1;
	c = c_values(item,c_values.columns(item)-1); 

	//Initialized with gamma value
	double eta = parameters(parameters.size() - 1);

	//Computing dot product
	for ( size_t i = 0; i < theta.size(); ++i )
		eta += parameters(i) * theta[i];

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

bayesian::~bayesian() {}

}

} /* namespace latentregpp */

