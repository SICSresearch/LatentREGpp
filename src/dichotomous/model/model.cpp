/*
 * model.cpp
 *
 *  Created on: 13/04/2016
 *      Author: Milder
 */

#include "../../dichotomous/model/model.h"

namespace latentregpp {

namespace dichotomous {

model::model(model_type type, int parameters) {
	this->type = type;
	this->parameters = parameters;
}

double model::P ( const optimizer_vector &theta, const optimizer_vector &parameters ) {
	std::vector<double> theta_temp(theta.size());
	for ( int i = 0; i < theta.size(); ++i )
		theta_temp[i] = theta(i);
	return P(theta_temp, parameters);
}

model::~model() {
	// TODO Auto-generated destructor stub
}

}

} /* namespace latentregpp */
