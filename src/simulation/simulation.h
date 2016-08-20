/*
 * simulation.h
 *
 *  Created on: 1/06/2016
 *      Author: Milder
 */

#ifndef SIMULATION_SIMULATION_H_
#define SIMULATION_SIMULATION_H_

#include <iostream>
#include "../util/matrix.h"
#include "../util/input.h"
#include "../util/constants.h"
#include "../util/quadraturepoints.h"
#include "../polytomous/estimation/estimation.h"
#include "../dicho-multi/estimation/estimation.h"
#include <fstream>
#include <cstdlib>
#include <ctime>
#include <cstdio>
#include <vector>
#include <set>
#include <omp.h>

namespace irtpp {

/**
 * Simulation class has a collection of helper functions
 * to simulate and testing dataset, model, package, precision, time
 * and other statistics from IRT.
 */
class simulation {
public:

	/**
	 * Helper function to simulate a specific range of the total number of iterations
	 * 		range -> [start, end], where 1 <= start <= iterations
	 * 									 1 <= end <= iterations
	 * 									 start <= end
	 *
	 * Should be called from.
	 * 		void simulate ( int model, int d, int iterations, std::string folder, std::string name, int interval );
	 *
	 * 		and it's called iterations/interval times
	 *
	 * @param model integer number with model, 1 = "1PL", 2 = "2PL", 3 = "3PL".
	 * @param d the dimension.
	 * @param start the start iteration.
	 * @param end the end iteration.
	 * @param folder string path where are the datasets.
	 * @param name relative path within folder with dataset(s) to estimate.
	 * @param dif epsilon convergence.
	 * @param dicho if true, the estimation is dichotomous.
	 * @param quadrature_technique it can be "Gaussian" or "Sobol".
	 * @param G number of points in Sobol technique.
	 * @param cluster a std vector cluster with integer template for multidimensional cases.
	 * @param custom_initial_values_filename BUILD if the dataset has custom initial values, else NONE.
	 */
	 void simulate ( int model, int d, int start, int end, std::string folder, std::string name, double, bool,
					std::string, int, std::vector<int> cluster, std::string );

	/**
	 * Simulates the number of iterations.
	 *
	 * And saves estimations each 'interval' iterations
	 *
	 *	Receives the model to be used and dimension of the problem,
	 *	number of iterations, the folder where the datasets will be loaded and estimation saved
	 *	and the prefix of files allocated in 'folder'
	 *	bool to indicate if data is dichotomous or not
	 *
	 *	Example:
	 *		Model: 2PL
	 *		Dimension: 2
	 *		Iterations: 100
	 *		Folder where data is: "datasets/dicho-multi-tests/escenario2"
	 *		Prefix of files in folder: "dicho-multi"
	 *		Convergence difference
	 *		Saves each 20 iterations.
	 *		Data is dichotomous
	 *
	 *		simulate(2, 2, 100, "datasets/dicho-multi-tests/escenario2", "dicho-multi", 20, 0.001, true);
	 *
	 *	Optional parameters:
	 *		quadrature_technique: [GAUSSIAN_QUADRATURE, SOBOL_QUADRATURE]
	 *		G: 					  Custom number of quadrature points
	 *		cluster: 			  Vector that contains the number of item for each dimension
	 *		custom_initial_values_filename:		  Custom initial values filename
	 *
	 * @param model integer number with model, 1 = "1PL", 2 = "2PL", 3 = "3PL".
	 * @param d the dimension.
	 * @param iterations the number of iterations.
	 * @param folder string path where are the datasets.
	 * @param name relative path within folder with dataset(s) to estimate.
	 * @param interval save the estimations each interval number.
	 * @param dif epsilon convergence.
	 * @param dicho if true, the estimation is dichotomous.
	 * @param quadrature_technique it can be "Gaussian" or "Sobol".
	 * @param G number of points in Sobol technique.
	 * @param cluster a std vector cluster with integer template for multidimensional cases.
	 * @param custom_initial_values_filename BUILD if the dataset has custom initial values, else NONE.
	 */
	void simulate ( int model, int d, int iterations, std::string folder,
				    std::string name, int interval, double, bool,
					std::string quadrature_technique = QMCEM,
					int G = DEFAULT_QMCEM_POINTS,
					std::vector<int> cluster = std::vector<int>(),
					std::string custom_initial_values_filename = NONE );

	/**
	 * Runs a single test.
	 * Receives the model to use
	 * 			the dimensions of the problem
	 * 			the filename of the dataset
	 * 			the convergence difference to use
	 * 			bool to indicate if data is dichotomous or not
	 *
	 * Optional parameters:
	 *		quadrature_technique: [GAUSSIAN_QUADRATURE, SOBOL_QUADRATURE]
	 *		G: 					  Custom number of quadrature points
	 *		cluster: 			  Vector that contains the number of item for each dimension
	 *		custom_initial_values_filename:		  Custom initial values filename
	 *
	 * This function calls dichotomous or polytomous according bool variable
	 *
	 * @param model integer number with model, 1 = "1PL", 2 = "2PL", 3 = "3PL".
	 * @param d the dimension.
	 * @param filename string path where is the dataset.
	 * @param dif epsilon convergence.
	 * @param dicho if true, the estimation is dichotomous.
	 * @param quadrature_technique it can be "Gaussian" or "Sobol".
	 * @param G number of points in Sobol technique.
	 * @param cluster a std vector cluster with integer template for multidimensional cases.
	 * @param custom_initial_values_filename BUILD if the dataset has custom initial values, else NONE.
	 */
	void run_single ( int, int, std::string, double, bool,
					  std::string quadrature_technique = GAUSSIAN, int G = DEFAULT_QMCEM_POINTS,
					  std::vector<int> cluster = std::vector<int>(),
					  std::string custom_initial_values_filename = NONE );

	/**
	 * Runs a single polytomous test
	 * Receives the model to use
	 * 			the dimensions of the problem
	 * 			the filename of the dataset
	 * 			the convergence difference to use
	 *
	 *	Optional parameters:
	 *		quadrature_technique: [GAUSSIAN_QUADRATURE, SOBOL_QUADRATURE]
	 *		G: 					  Custom number of quadrature points
	 *		cluster: 			  Vector that contains the number of item for each dimension
	 *		custom_initial_values_filename:		  Custom initial values filename
	 *
	 * @param model integer number with model, 1 = "1PL", 2 = "2PL", 3 = "3PL".
	 * @param d the dimension.
	 * @param filename string path where is the dataset.
	 * @param dif epsilon convergence.
	 * @param quadrature_technique it can be "Gaussian" or "Sobol".
	 * @param G number of points in Sobol technique.
	 * @param cluster a std vector cluster with integer template for multidimensional cases.
	 * @param custom_initial_values_filename BUILD if the dataset has custom initial values, else NONE.
	 */
	void run_single_polytomous ( int, int, std::string, double,
								 std::string quadrature_technique = GAUSSIAN, int G = DEFAULT_QMCEM_POINTS,
								 std::vector<int> cluster = std::vector<int>(),
								 std::string custom_initial_values_filename = NONE );

	/**
	 * Runs a single polytomous test
	 * Receives the model to use
	 * 			the dimensions of the problem
	 * 			the filename of the dataset
	 * 			the convergence difference to use
	 *
	 *	Optional parameters:
	 *		quadrature_technique: [GAUSSIAN_QUADRATURE, SOBOL_QUADRATURE]
	 *		G: 					  Custom number of quadrature points
	 *		cluster: 			  Vector that contains the number of item for each dimension
	 *		custom_initial_values_filename:		  Custom initial values filename
	 *
	 * @param model integer number with model, 1 = "1PL", 2 = "2PL", 3 = "3PL".
	 * @param d the dimension.
	 * @param filename string path where is the dataset.
	 * @param dif epsilon convergence.
	 * @param quadrature_technique it can be "Gaussian" or "Sobol".
	 * @param G number of points in Sobol technique.
	 * @param cluster a std vector cluster with integer template for multidimensional cases.
	 * @param custom_initial_values_filename BUILD if the dataset has custom initial values, else NONE.
	 */
	void run_single_dichotomous ( int, int, std::string, double,
								  std::string quadrature_technique = GAUSSIAN, int G = DEFAULT_QMCEM_POINTS,
								  std::vector<int> cluster = std::vector<int>(),
								  std::string custom_initial_values_filename = NONE );

	/**
	 * Default constructor for simulation class. Is not working, don't use.
	 */
	simulation();

	/**
	 * Default destructor for simulation class. Is not working, don't use.
	 */
	virtual ~simulation();
};

} /* namespace irtpp */

#endif /* SIMULATION_SIMULATION_H_ */
