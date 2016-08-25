/**
 * estep.h
 *
 *  Created on: 13/04/2016
 *      Author: Milder
 */

#ifndef POLYTOMOUS_ESTIMATION_ESTEP_H_
#define POLYTOMOUS_ESTIMATION_ESTEP_H_

#include "../model/model.h"

#include "../../util/matrix.h"
#include "../type/estimationdata.h"

#include "../../test/test.h"
#include <ctime>

#include <iostream>

#include <omp.h>

namespace latentregpp {

namespace polytomous {


/**
 * Estep of the EMAlgorithm.
 *
 * Receives an estimation_data reference that MUST bring all the
 * data needed to run the Estep.
 * @param data a instance of estimation_data with all necessary data to do EMAlgorithm.
 */
void Estep(estimation_data&, int current);

} /* namespace latentregpp */

}

#endif /* ESTIMATION_ESTEP_H_ */
