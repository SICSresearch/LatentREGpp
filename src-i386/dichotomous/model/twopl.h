/*
 * twopl.h
 *
 *  Created on: Sep 20, 2016
 *      Author: milder
 */

#ifndef DICHOTOMOUS_MODEL_TWOPL_H_
#define DICHOTOMOUS_MODEL_TWOPL_H_

#include "model.h"

namespace latentregpp {

namespace dichotomous {

class twopl: public dichotomous::model {
public:
	twopl();
	virtual ~twopl();
	double P(std::vector<double>&, const optimizer_vector&);
};

}

} /* namespace latentregpp */

#endif /* DICHOTOMOUS_MODEL_TWOPL_H_ */
