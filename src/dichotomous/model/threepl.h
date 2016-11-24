/*
 * threepl.h
 *
 *  Created on: Sep 20, 2016
 *      Author: milder
 */

#ifndef DICHOTOMOUS_MODEL_THREEPL_H_
#define DICHOTOMOUS_MODEL_THREEPL_H_

#include "model.h"

namespace latentregpp {

namespace dichotomous {

class threepl: public dichotomous::model {
public:
	threepl();
	virtual ~threepl();
	double P(std::vector<double>&, const optimizer_vector&, int);
};

}

} /* namespace latentregpp */

#endif /* DICHOTOMOUS_MODEL_THREEPL_H_ */
