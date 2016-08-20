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

#include "../../util/matrix.h"
#include "../../util/input.h"
#include "../../util/quadraturepoints.h"
#include "../../util/initial_values.h"
#include "../../util/constants.h"
#include "../../util/ramsay.h"

#include "../type/estimationdata.h"

#include <map>
#include <cmath>

namespace irtpp {

namespace polytomous {

/**
 * Class to set up and run the estimation process
 * in polytomous models
 *
 * The main method is EMAlgorithm
 */
class estimation {
	private:

		short iterations; /**< Counts the actual number of iterations*/
		double convergence_difference; /**< Epsilon to stop the EMAlgorithm*/
		std::string custom_initial_values_filename; /**< path with custom initial values*/

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
		 * Receives:
		 * 	1. A specific model to use -> [1, 3] that means 1PL, 2PL, or 3PL
		 *  2. A matrix containing the answers of examinees
		 *  3. The dimension of the problem
		 *	4. The epsilon (convergences difference) that the algoritm will use
		 *		as a criterion to stop
		 *
		 *	Optional parameters:
		 *		quadrature_technique: [GAUSSIAN_QUADRATURE, SOBOL_QUADRATURE]
		 *		quadrature_points: 	  Custom number of quadrature points
		 *		cluster: 			  Vector that contains the number of item for each dimension
		 *		custom_initial_values_filename:		  Custom initial values filename
		 *
		 *  Then it sets up all the data needed to start the estimation process
		 *
		 * @param dataset a matrix data type char template with data to estimate parameters.
		 * @param d the dimension.
		 * @param themodel model to use 1PL, 2PL or 3PL.
		 * @param convergence_difference epsilon convergence difference.
		 * @param quadrature_technique string. it can be Gaussian or Sobol.
		 * @param quadrature_technique if Sobol. Number of points to use.
		 * @param cluster a std vector integer template with number of items for each dimension.
		 * @param custom_initial_values_filename string with path for custom initial_values. Default is none.
		 */
		estimation(matrix<char>&, unsigned int, int themodel = 2, double convergence_difference = 0.001,
								std::vector<int> cluster = EMPTY_INTEGER_VECTOR,
								std::string quadrature_technique = GAUSSIAN,
								int quadrature_points = DEFAULT_QMCEM_POINTS,
								std::vector<int> individuals_weights = EMPTY_INTEGER_VECTOR,
								std::string custom_initial_values_filename = NONE );

		/**
		 * Destructor for estimation class.
		 */
		virtual ~estimation();

		/**
		 * Finds the initial values for every parameter of the items to start the estimation.
		 * it is called if custom_initial_values_filename is none by default.
		 * @see custom_initial_values_filename
		 */
		void initial_values();

		/**
		 * Loads the initial values for every parameter of the items to start the estimation
		 * from file. It is call if custom_initial_values_filename has a value different to none.
		 * @param filename string with path for initial values.
		 * @see custom_initial_values_filename
		 */
		void load_initial_values(std::string);

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
		void EMAlgortihm();

		/*
		 * EAP
		 * */
		void EAP(bool);

		/*
		 * MAP
		 * */
		void MAP(bool);

		/**
		 * Prints the item parameters values of the estimated parameters.
		 */
		void print_item_parameters();

		/**
		 * Prints the item parameters to a specific file, including time elapsed in ms.
		 * @param filename the filename of the where parameters will be saved
		 * @param elapsed a double with time elapsed in EM estimation.
		 */
		void print_item_parameters(std::string, double elapsed = -1.0);

		/**
		 * Prints the results to a specific output stream, including time elapsed in ms.
		 * @param fout the ofstream object from std library.
		 * @param elapsed a double with time elpased in EM estimation.
		 */
		void print_item_parameters(std::ofstream &, double elapsed = -1.0);

		/**
		 * Prints the latent traits.
		 */
		void print_latent_traits();

		/**
		 * Prints the latent traits to a specific file
		 * @param filename the filename of the where latent traits will be saved
		 */
		void print_latent_traits(std::string);

		/**
		 * Prints the latent traits to a specific output stream
		 * @param fout the ofstream object from std library.
		 */
		void print_latent_traits(std::ofstream &);

		/**
		 * @return log_likehood of the estimation
		 * */
		double log_likelihood();

		/**
		 * @return number of iterations
		 * */
		short get_iterations();
};

}

} /* namespace irtpp */

#endif /* ESTIMATION_ESTIMATION_H_ */
