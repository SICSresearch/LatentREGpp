/*
 * onepl.h
 *
 *  Created on: Sep 20, 2016
 *      Author: milder
 */

#ifndef DICHOTOMOUS_MODEL_ONEPL_H_
#define DICHOTOMOUS_MODEL_ONEPL_H_

#include "model.h"

namespace latentregpp {

namespace dichotomous {

class onepl: public dichotomous::model {
public:
	onepl();
	virtual ~onepl();
	double P(std::vector<double>&, const optimizer_vector&);
};

}

} /* namespace latentregpp */

#endif /* DICHOTOMOUS_MODEL_ONEPL_H_ */
