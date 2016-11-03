/*
 * bayesian.h
 *
 *  Created on: 26/10/2016
 *      Author: jhonatan
 */

#ifndef DICHOTOMOUS_MODEL_BAYESIAN_H_
#define DICHOTOMOUS_MODEL_BAYESIAN_H_

#include "model.h"

namespace latentregpp {

namespace dichotomous {

class bayesian: public dichotomous::model {
public:
	//bayesian();
	bayesian(matrix<double>&);
	virtual ~bayesian();
	double P(std::vector<double>&, const optimizer_vector&, int);

private:
	matrix<double> c_values; 
};

}

} /* namespace latentregpp */

#endif /* DICHOTOMOUS_MODEL_BAYESIAN_H_ */
