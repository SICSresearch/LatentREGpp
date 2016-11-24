/**
 * estep.h
 *
 *  Created on: 13/04/2016
 *      Author: Milder
 */

#ifndef dichotomous_ESTIMATION_ESTEP_H_
#define dichotomous_ESTIMATION_ESTEP_H_

#include "../../util/matrix.h"
#include "../../test/test.h"
#include <ctime>

#include <iostream>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "../../dichotomous/model/model.h"
#include "../../dichotomous/type/estimationdata.h"

namespace latentregpp {

namespace dichotomous {

/**
 * Estep of the EMAlgorithm.
 *
 * Receives an estimation_data reference that MUST bring all the
 * data needed to run the Estep
 * @param data a instance of estimation_data with all necessary data to do EMAlgorithm.
 * @param current the current zeta estimation (Ramsay and Squarem accelerate)
 */
void Estep(estimation_data&, int current);

} /* namespace latentregpp */

}

#endif /* ESTIMATION_ESTEP_H_ */
