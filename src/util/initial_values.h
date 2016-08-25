/**
 * initial_values.h
 *
 *  Created on: 13/06/2016
 *      Author: Milder
 */

#ifndef UTIL_INITIAL_VALUES_H_
#define UTIL_INITIAL_VALUES_H_

#include <vector>
#include <exception>
#include "matrix.h"

namespace latentregpp {
	/**
	 * Computing initial values according to
	 * Andrade, Tavares & Valle (2000),
	 *
	 */


	/**
	 * Computes initial values for Unidimensional - dichotomous case.
	 * @param data a matrix data type with char template with initial data by reference.
	 * @param alpha a std vector with double template where alpha's initial values are saved, pass by reference.
	 * @param gamma a std vector with double template where gamma's initial values are saved, pass by reference.
	 */
	void find_initial_values ( matrix<char> &data, std::vector<double> &alpha, std::vector<double> &gamma );

	/**
	 * Computes biserial correlation.
	 * Needed to find a(discrimation) and b(difficulty) initial values.
	 *
	 * PARAMETERS:
	 * 		IN: scores, Y
	 * 		OUT: p  ---> biserial correlation of each item
	 *
	 * p is Ro greek letter in initial values document.
	 *
	 * @param Y matrix with char template, pass by reference.
	 * @param p a std vector with double template where will be saved biserial correlation, pass by reference.
	 */
	void compute_biserial_correlation(matrix<char> &Y, std::vector<double> &p);

	/**
	 * Compute alpha's initial values.
	 * p is Ro greek letter in initial values document.
	 * @param p std vector with double template as biserial correlation, pass by reference.
	 * @param a std vector with double template where alphas will be saved.
	 */
	void compute_alphas(std::vector<double> &p, std::vector<double> &a);

	/**
	 * Compute gamma's initial values.
	 * @param Y matrix with char template, pass by reference.
	 * @param a std vector with double template as alpha values, pass by reference.
	 * @param p std vector with double template as biserial correlation, pass by reference.
	 * @param d std vector with double template where gammas will be saved, pass by reference.
	 */
	void compute_gammas(matrix<char> &Y, std::vector<double> &a, std::vector<double> &p, std::vector<double> &d);


	/**
	 * Helper function to find mean of an array.
	 * @param v a std vector with double template to calculate mean.
	 * @return the mean of v vector
	 */
	double mean ( std::vector<double> & );

	/**
	 * Helper function to find standard deviation of an array.
	 * @param v a std vector with double template to calculate standard deviation.
	 * @return the standard deviation of v vector
	 */
	inline double sd ( std::vector<double> &v );
} //latentregpp

#endif /* UTIL_INITIAL_VALUES_H_ */
