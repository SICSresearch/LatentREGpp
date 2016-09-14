/**
 * mstep.h
 *
 *  Created on: 13/04/2016
 *      Author: Milder
 */

#ifndef dichotomous_ESTIMATION_MSTEP_H_
#define dichotomous_ESTIMATION_MSTEP_H_

#include "../../util/matrix.h"
#include <iostream>
#include <vector>
#include <cmath>
#include <set>

#include "../../dichotomous/model/model.h"
#include "../../dichotomous/type/estimationdata.h"

namespace latentregpp {

namespace dichotomous {

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
double Mstep(estimation_data&, int);

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
  * @param current_zeta_obj the current zeta
  */
  Qi (int, estimation_data*, int);
  
  /**
   * Overload parenthesis operator to evaluate
   * the Log likelihood function.
   */
  double operator() (const optimizer_vector&) const;
  
private:
  int i; /**< The current item*/
  estimation_data *data; /**< estimation_data pointer*/
  int current_zeta_obj;
};

}

} /* namespace latentregpp */
  
#endif /* ESTIMATION_MSTEP_H_ */
  