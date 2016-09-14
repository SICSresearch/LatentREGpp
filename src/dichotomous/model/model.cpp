/*
 * model.cpp
 *
 *  Created on: 13/04/2016
 *      Author: Milder
 */

#include "../../dichotomous/model/model.h"

namespace latentregpp {

namespace dichotomous {

model::model() : model(2) {
}

model::model(int parameters) {
	this->parameters = parameters;
}

double model::P ( const optimizer_vector &theta, const optimizer_vector &parameters ) {
	std::vector<double> theta_temp(theta.size());
	for ( int i = 0; i < theta.size(); ++i )
		theta_temp[i] = theta(i);
	return P(theta_temp, parameters);
}

double model::P(std::vector<double> &theta, const optimizer_vector &parameters) {
	if ( this->parameters == 1 ) {
		/**
		 * 1PL Approach
		 *
		 * */

		//Initialized with gamma_k value
		double eta = parameters(0);

		//Computing dot product
		for ( size_t i = 0; i < theta.size(); ++i )
			eta += 1 * theta[i]; //no alpha in this model

		double P = 1.0 / (1.0 + std::exp(-eta));
		P = std::max(P, LOWER_BOUND_);
		P = std::min(P, UPPER_BOUND_);
		return P;
	}
	if ( this->parameters == 2 ) {
		/**
		 * 2PL Approach
		 * */

		//Initialized with gamma value
		double eta = parameters(parameters.size() - 1);

		//Computing dot product
		for ( size_t i = 0; i < theta.size(); ++i )
			eta += parameters(i) * theta[i];

		double P = 1.0 / (1.0 + std::exp(-eta));
		P = std::max(P, LOWER_BOUND_);
		P = std::min(P, UPPER_BOUND_);
		return P;
	}
	/**
	 * 3PL Approach
	 * */
	double gamma_parameter = parameters(parameters.size() - 1);
	
	//uncommented line below for reparameter a c value [0,1] from gamma in R
	double c = 1.0 / (1.0 + exp(-gamma_parameter));

	//Initialized with gamma value
	double eta = parameters(parameters.size() - 2);

	//Computing dot product
	for ( size_t i = 0; i < theta.size(); ++i )
		eta += parameters(i) * theta[i];

	/**three differente formulas**/
	
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

model::~model() {
	// TODO Auto-generated destructor stub
}

}

} /* namespace latentregpp */
