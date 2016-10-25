/**
 * quadrature_points.h
 *
 *  Created on: 22/04/2016
 *      Author: Milder
 */

#ifndef SRC_UTIL_QUADRATUREPOINTS_H_
#define SRC_UTIL_QUADRATUREPOINTS_H_

#include <cstring>
#include <bitset>
#include <iostream>
#include <fstream>
#include <sstream>
#include "matrix.h"
#include "directories.h"
#include <vector>
#include <cmath>

namespace latentregpp {
	/**
	 * Functions to compute and save quadrature points and weights
	 * based on gauss.quad() function from statmod library in R
	 * */

	/**
	 * Function to compute and save quadrature points.
	 * will be saved within data folder with .data extension
	 * @param G an integer with number of points
	 * @param d an integer with the dimension
	 */
	void compute_and_save_quadrature_points(int, int);

	/**
	 * Function to compute and save weights.
	 * will be saved within data folder with .data extension
	 * @param G an integer with number of points
	 * @param d an integer with the dimension
	 */
	void compute_and_save_weights(int, int);


	/**
	 * Functions to load quadrature points previously computed
	 * by functions above
	 * @param d the dimension for load those points
	 * @return a matrix data type double template with quad points loaded
	 * */
	matrix<double> load_quadrature_points(int);

	/**
	 * Function to load weights previously computed
	 * by functions above
	 * @param d the dimension for load those points
	 * @return a std vector double template with weights loaded
	 */
	std::vector<double> load_weights(int);
}

#endif /* SRC_UTIL_QUADRATUREPOINTS_H_ */
