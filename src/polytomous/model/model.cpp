/*
 * model.cpp
 *
 *  Created on: 13/04/2016
 *      Author: Milder
 */

#include "model.h"

namespace latentregpp {

namespace polytomous {

model::model() {}

model::model(model_type type, int parameters, int d, std::vector<int> *categories_item) {
	this->type = type;
	this->parameters = parameters;
	this->d = d;
	this->categories_item = categories_item;
}

model::~model() {}

double model::Pik(const optimizer_vector &theta, const optimizer_vector &parameters, int i, int k) {
	std::vector<double> theta_temp(theta.size());
	for ( int h = 0; h < theta.size(); ++h )
		theta_temp[h] = theta(h);
	return Pik(theta_temp, parameters, i, k);
}

}

} /* namespace latentregpp */
