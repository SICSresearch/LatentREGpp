/*
 * estimationdata.h
 *
 *  Created on: 20/06/2016
 *      Author: Milder
 */

#ifndef POLYTOMOUS_UTIL_ESTIMATIONDATA_H_
#define POLYTOMOUS_UTIL_ESTIMATIONDATA_H_
#include <vector>
#include <set>
#include "../model/model.h"
#include "../../util/matrix.h"
#include "../../util/constants.h"

#include <dlib/optimization.h>

namespace latentregpp {

namespace polytomous {

/**
 * estimation_data class contains all the information needed to execute the EM estimation in
 * polytomous model
 */
class estimation_data {
public:
	matrix<char> *dataset; /**< Matrix of answers*/
	int d; /**< Dimension*/
	matrix<char> Y; /**< Matrix of response patterns*/
	std::vector<int> nl; /**< Frequencies of each response pattern*/
	int N; /**< Number of examines*/
	int s; /**< Number of response patterns*/
	int p; /**< Number of items*/
	int G; /**< Number of quadrature points*/
	std::vector<int> categories_item; /**< Number of categories for each item*/
	matrix<double> theta; /**< Latent traits vectors*/
	std::vector<double> w; /**< Weights*/
	std::vector<matrix<double> > r; /**< Matrix r*/
	std::vector<matrix<double> > P; /**< Probability matrix P, P_gik means the probability that an individual has selected the category k to item i and belongs to group g*/
	matrix<double> pi; /**< Matrix pi*/
	std::set<int> pinned_items; /**< Pinned items (won't be estimated)*/
	model *m; /**< Model to use*/
	double loglikelihood; /**< Loglikelihood */
	std::vector<optimizer_vector> zeta[ACCELERATION_PERIOD]; /**< Vector of zeta item parameters*/
	std::map<std::vector<char>, std::vector<int> > patterns; /**< Patterns and their individuals*/
	std::vector<optimizer_vector> latent_traits; /**< Latent traits */

	/**
	 * Constructor for estimation_data class.
	 * @param d the dimension.
	 */
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
