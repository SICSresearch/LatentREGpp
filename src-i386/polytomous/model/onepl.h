/*
 * onepl.h
 *
 *  Created on: Sep 20, 2016
 *      Author: milder
 */

#ifndef POLYTOMOUS_MODEL_ONEPL_H_
#define POLYTOMOUS_MODEL_ONEPL_H_

#include "model.h"

namespace latentregpp {

namespace polytomous {

class onepl: public polytomous::model {
public:
	onepl();
	virtual ~onepl();
	onepl(int, std::vector<int>*);

	/**
	 * Probability in dichotomous case according to equation 82 in IRT_engineers document.
	 * @param theta a std vector double template with theta values.
	 * @param parameters a typedef optimizer_vector to extract eta values.
	 * @param i number of items.
	 * @param k number of categories.
	 * @return the probability in dichotomous case according select model. it can be 1PL or 2PL
	 */
	double Pstar_ik(std::vector<double>&, const optimizer_vector&, int i, int k);

	double Pik(std::vector<double>&, const optimizer_vector&, int i, int k);
};

}

} /* namespace latentregpp */

#endif /* POLYTOMOUS_MODEL_ONEPL_H_ */
