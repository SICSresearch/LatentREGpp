/**
 * mstep.h
 *
 *  Created on: 13/04/2016
 *      Author: Milder
 */

#ifndef DICHOMULTI_ESTIMATION_MSTEP_H_
#define DICHOMULTI_ESTIMATION_MSTEP_H_

#include "../../util/matrix.h"
#include <iostream>
#include <vector>
#include <cmath>
#include <set>

#include "../../dicho-multi/model/model.h"
#include "../../dicho-multi/type/estimationdata.h"

namespace irtpp {

namespace dichomulti {

/**
 * M step of the EM Algorithm.
 *
 * Receives an estimation_data reference that MUST bring all the
 * data needed to run the Mstep.
 *
 * @param data a instance of estimation_data with all necessary data to do EMAlgorithm.
 * @param current_zeta the current zeta estimation (Ramsay and Squarem accelerate).
 * @return the max difference between current zeta optimized and previous zeta
 */
double Mstep(estimation_data&, int current_zeta);

/**
 * Class with Log likelihood Function to maximize
 */
class Qi {
public:
	/**
	 * Constructor that receives the number of the current item (i)
	 * and the estimation_data pointer
	 * @param i the current item
	 * @param d estimation_data pointer
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

}

} /* namespace irtpp */

#endif /* ESTIMATION_MSTEP_H_ */
