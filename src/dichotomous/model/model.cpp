/*
 * model.cpp
 *
 *  Created on: 13/04/2016
 *      Author: Milder
 */

#include "../../dichotomous/model/model.h"

namespace latentregpp {

namespace dichotomous {

model::model() {}

model::model(model_type type, int parameters) {
	this->type = type;
	this->parameters = parameters;
}

double model::P ( const optimizer_vector &theta, const optimizer_vector &parameters, int i ) {
	std::vector<double> theta_temp(theta.size());
	for ( int j = 0; j < theta.size(); ++j )
		theta_temp[j] = theta(j);
	return P(theta_temp, parameters, i);
}

model::~model() {}

}

} /* namespace latentregpp */
