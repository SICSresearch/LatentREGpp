/**
 * model.h
 *
 *  Created on: 13/04/2016
 *      Author: Milder
 */

#ifndef POLYTOMOUS_MODEL_MODEL_H_
#define POLYTOMOUS_MODEL_MODEL_H_

#include <vector>
#include <cmath>
#include <cassert>

#include "../../util/constants.h"

namespace latentregpp {

namespace polytomous {

/**
 * Model class
 * It represents what is the model approach to use
 * Might be 1PL or 2PL.
 */
class model {

public:
	model_type type; /**< model type */
	int parameters; /**< Number of parameters of the model*/
	std::vector<int> *categories_item; /**< Number of categories for each item*/
	int d; /**< the dimension*/

	/**
	 * Default constructor for model class. It's empty, don't use.
	 */
	model();

	/**
	 * Constructor. This receives 1, 2 or 3. Depending on the model to use, dimension and
	 * number of categories for each item.
	 * @param type an model_type that represents the type of the model
	 * @param parameters an integer with number of parameters. Decide the model 1PL, 2PL or 3PL.
	 * @param d the dimension.
	 * @param categories_item a std vector integer template with number of categories for each item.
	 */
	model(model_type, int, int, std::vector<int>*);

	/**
	 * This method computes the probability that a response pattern U_l has the category k to item
	 * i, given that its latent trait vector theta, and the item parameters. According to equation 89
	 * in IRT_engineers document.
	 * @param theta a std vector double template with theta values.
	 * @param parameters a typedef optimizer_vector to extract eta values.
	 * @param i number of items.
	 * @param k number of categories.
	 * @return the probability for polytomous.
	 */
	double Pik(const optimizer_vector&, const optimizer_vector&, int i, int k);


	/**
	 * This method computes the probability that a response pattern U_l has the category k to item
	 * i, given that its latent trait vector theta, and the item parameters. According to equation 89
	 * in IRT_engineers document.
	 * @param theta a std vector double template with theta values.
	 * @param parameters a typedef optimizer_vector to extract eta values.
	 * @param i number of items.
	 * @param k number of categories.
	 * @return the probability for polytomous.
	 */
	virtual double Pik(std::vector<double>&, const optimizer_vector&, int i, int k) = 0;

	/**
	 * Destructor for model class.
	 */
	virtual ~model();
};

}

} /* namespace latentregpp */

#endif /* MODEL_MODEL_H_ */
