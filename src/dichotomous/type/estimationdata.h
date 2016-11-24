/*
 * estimationdata.h
 *
 *  Created on: 20/06/2016
 *      Author: Milder
 */

#ifndef dichotomous_UTIL_ESTIMATIONDATA_H_
#define dichotomous_UTIL_ESTIMATIONDATA_H_
#include <vector>
#include <set>
#include "../../util/matrix.h"
#include "../../util/constants.h"

#include <dlib/optimization.h>

#include <algorithm>

#include "../../dichotomous/model/model.h"

namespace latentregpp {

namespace dichotomous {

/**
 * estimation_data class contains all the information needed to execute the EM estimation in
 * dichotomous model
 */
class estimation_data {
public:
	matrix<char> *dataset; /**< Matrix of answers*/
	matrix<int> correct; /**< Matrix that contains what items have been answered correctly for each response pattern*/
	int d; /**< Dimension*/
	matrix<char> Y; /**< Matrix of response patterns*/
	std::vector<double> nl; /**< Frequencies of each response pattern*/
	int N; /**< Number of examines*/
	int s; /**< Number of response patterns*/
	int p; /**< Number of items*/
	int G; /**< Number of quadrature points*/
	matrix<double> theta; /**< Latent traits vectors*/
	std::vector<double> w; /**< Weights*/
	matrix<double> r; /**< Matrix r*/
	matrix<double> P; /**< Probability matrix P, P_gi means the probability that an individual has selected the correct answer*/
	matrix<double> pi; /**< Matrix pi*/
	std::vector<double> f; /**< Vector f (Number of individuals in group g)*/
	std::set<int> pinned_items; /**< Pinned items (won't be estimated)*/
	model *m; /**< Model to use*/
	double loglikelihood; /**< Loglikelihood */
	std::vector<optimizer_vector> zeta[ACCELERATION_PERIOD]; /**< Vector of zeta item parameters*/
	std::map<std::vector<char>, std::vector<int> > patterns; /**< Patterns and their individuals*/
	std::vector<optimizer_vector> latent_traits; /**< Latent traits */
  	matrix<double> initial_values; /**< A vector with initial values for data */
  	bool noguessing; /**< Flag for bayesian models, decision for remove or allow c parameter */

	estimation_data(int);

	/**
	 * Default constructor for estimation_data. Is not working, don't use.
	 */
	estimation_data();

	/**
	 * Default destructor for estimation_data. Is not working, don't use.
	 */
	virtual ~estimation_data();
};

}

} /* namespace latentregpp */

#endif /* UTIL_ESTIMATIONDATA_H_ */
