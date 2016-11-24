/**
 * estimation.h
 *
 *  Created on: 13/04/2016
 *      Author: Milder
 */

#ifndef POLYTOMOUS_ESTIMATION_ESTIMATION_H_
#define POLYTOMOUS_ESTIMATION_ESTIMATION_H_

#include "estep.h"
#include "mstep.h"

#include "../model/model.h"
#include "../model/onepl.h"
#include "../model/twopl.h"

#include "../../util/matrix.h"
#include "../../util/input.h"
#include "../../util/quadraturepoints.h"
#include "../../util/initial_values.h"
#include "../../util/constants.h"
#include "../../util/ramsay.h"

#include "../type/estimationdata.h"

#include <map>
#include <cmath>
#include <cstdio>

namespace latentregpp {

namespace polytomous {

/**
 * Class to set up and run the estimation process
 * in polytomous models
 *
 * The main method is EMAlgorithm
 */
class estimation {
	private:

		unsigned int iterations; /**< Counts the actual number of iterations*/
		double convergence_difference; /**< Epsilon to stop the EMAlgorithm*/

		/**
		 * Class to maximize posterior function to estimate latent traits
		 * in MAP approach
		 * */
		class posterior {
		public:
			/**
			 * Constructor that receives the number of the current item (i)
			 * and the estimation_data pointer
			 * @param l the current latent trait
			 * @param current_zeta the current zeta estimation
			 * @param d estimation_data pointer
			 */
			posterior (int, int, estimation_data*);

			/**
			 * Overload parenthesis operator to evaluate the function
			 */
			double operator() (const optimizer_vector&) const;
		private:
			int l; /**< The current latent trait*/
			int current_zeta; /**< The current zeta estimation (Ramsay and Squarem accelerate).*/
			estimation_data *data; /**< estimation_data pointer*/
		};

		/**
		 * Converts latent traits arranged by pattern to arranged by individuals
		 * */
		void latent_traits_by_individuals();

	public:

		estimation_data data; /**< Saves all data needed in the estimation process*/

		/**
		 *
		 * @param dataset a matrix data type char template with data to estimate parameters.
		 * @param d the dimension.
		 * @param themodel model to use 1PL, 2PL or 3PL.
		 * @param convergence_difference epsilon convergence difference.
		 * @param theta quadrature points
		 * @param weights quadrature points' weights
		 * @param pinned_items a std vector integer template with number of items for each dimension.
		 * @param initial_values matrix containing initial values
		 */
		estimation(matrix<char>&, unsigned int, int themodel = model_type::twopl,
						double convergence_difference = DEFAULT_EM_DELTA_STOP,
						matrix<double> theta = EMPTY_REAL_MATRIX,
					    std::vector<double> weights = EMPTY_REAL_VECTOR,
					    std::vector<double> individual_weights = EMPTY_REAL_VECTOR,
					    std::vector<int> pinned_items = EMPTY_INTEGER_VECTOR,
					    matrix<double> initial_values = EMPTY_REAL_MATRIX );

		/**
		 * Destructor for estimation class.
		 */
		virtual ~estimation();

		/**
		 * Finds the initial values for every parameter of the items to start the estimation.
		 * it is called if custom_initial_values_filename is none by default.
		 * @see custom_initial_values_filename
		 */
		void compute_initial_values();

		/**
		 * Loads the initial values for every parameter of the items to start the estimation
		 * @param mt matrix containing initial values
		 */
		void load_multi_initial_values(matrix<double> &mt);

		/**
		 * Sobol quadrature, receives the number of points to use.
		 * @param g the number of points.
		 */
		void sobol_quadrature (int);

		/**
		 * gaussian_quadrature. By default max points to use is 40.
		 */
		void gaussian_quadrature ();

		/**
		 * Builds all necessary matrixes for the estimation process
		 * */
		void build_matrixes();

		/**
		 * Runs the EMAlgorithm to find out the parameters.
		 */
		void EMAlgorithm(bool);

		/*
		 * EAP
		 * */
		void EAP(bool);

		/*
		 * MAP
		 * */
		void MAP(bool);

		/**
		 * @return log_likehood of the estimation
		 * */
		double log_likelihood();

		/**
		 * @return number of iterations
		 * */
		unsigned int get_iterations();

		void print_item_parameters();
};

}

} /* namespace latentregpp */

#endif /* ESTIMATION_ESTIMATION_H_ */
