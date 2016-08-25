/*
 * model.h
 *
 *  Created on: 13/04/2016
 *      Author: Milder
 */

#ifndef dichotomous_MODEL_MODEL_H_
#define dichotomous_MODEL_MODEL_H_

#include <vector>
#include <cmath>
#include <cassert>

#include "../../util/constants.h"

namespace latentregpp {

namespace dichotomous {

/**
 * Model class
 * It represents what is the model approach to use
 * Might be 1PL, 2PL, 3PL.
 */
class model {

public:
	int parameters; /**< Number of parameters of the model*/

	/**
	 * Default constructor for model class. It's empty, don't use.
	 */
	model();

	/**
	 * Constructor that receives 1, 2 or 3. Depending on the model to use.
	 * @param parameters number of parameters according model to use.
	 */
	model(int);

	/**
	 * Destructor for model class.
	 */
	virtual ~model();

	/**
	 * Function to calculate the probability according model to use and
	 * equation 13 from IRT_engineers document.
	 * @param theta a optimizer_vector with theta values.
	 * @param parameters a optimizer_vector to extract eta values.
	 * @return the probability given the model. It can be 1PL, 2PL or 3PL
	 */
	double P(const optimizer_vector &, const optimizer_vector &);

	/**
	 * Function to calculate the probability according model to use and
	 * equation 13 from IRT_engineers document.
	 * @param theta a std vector double template with theta values.
	 * @param parameters a optimizer_vector to extract eta values.
	 * @return the probability given the model. It can be 1PL, 2PL or 3PL
	 */
	double P(std::vector<double>&, const optimizer_vector&);
};

}

} /* namespace latentregpp */

#endif /* MODEL_MODEL_H_ */
