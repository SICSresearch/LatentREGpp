/*
 * mstep.h
 *
 *  Created on: 13/04/2016
 *      Author: Milder
 */

#ifndef POLYTOMOUS_ESTIMATION_MSTEP_H_
#define POLYTOMOUS_ESTIMATION_MSTEP_H_

#include "../model/model.h"

#include "../../util/matrix.h"

#include "../type/estimationdata.h"

#include <iostream>
#include <vector>
#include <cmath>
#include <set>

//including optimization files from dlib library
#include <dlib/optimization.h>

#include <omp.h>

namespace latentregpp {

namespace polytomous {

/**
 * M step of the EM Algorithm
 *
 * Receives an estimation_data reference that MUST bring all the
 * data needed to run the Mstep.
 * @param data a instance of estimation_data with all necessary data to do EMAlgorithm.
 * @return the max difference between current zeta optimized and previous zeta
 */
double Mstep(estimation_data&, int current);

/**
 * Class with Log likelihood Function to maximize
 */
class Qi {
public:
	/**
	 * Constructor that receives the number of the current item (i)
	 * and the estimation_data pointer.
	 * @param i the current item.
	 * @param d estimation_data pointer.
	 */
    Qi (int, estimation_data*);

    /**
	 * Overload parenthesis operator to evaluate
	 * the Log likelihood function.
	 */
    double operator() (const optimizer_vector&) const;
private:
    int i; /**< The current item*/
    estimation_data *data; /**< estimation_data pointer*/
};

/**
 * Class with derivatives for Log likelihood function
 */
class Qi_derivative {
public:
	/**
	 * constructor that receives the number of the current item (i)
	 * and the estimation_data pointer.
	 * @param i the current item.
	 * @param d estimation_data pointer.
	 */
	Qi_derivative (int, estimation_data*);

	/**
	 * Overload parenthesis operator to calculate lambda derivative
	 * and k derivatives for Log likelihood function.
	 */
	const optimizer_vector operator() (const optimizer_vector&) const;
private:
    int i; /**< The current item*/
    estimation_data *data; /**< estimation_data pointer*/
};

}

} /* namespace latentregpp */

#endif /* ESTIMATION_MSTEP_H_ */
