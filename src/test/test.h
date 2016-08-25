/*
 * pimatrixtest.h
 *
 *  Created on: 25/04/2016
 *      Author: Milder
 */

#ifndef TEST_TEST_H_
#define TEST_TEST_H_

#include <vector>
#include "../util/matrix.h"

namespace latentregpp {

	/**
	 * Tests the sum if each columns of pi matrix is equals to 1.
	 * @param pi a matrix data type with double template.
	 * @return true if pi matrix is correct
	 */
	bool test_pi ( matrix<double> &pi );

	/**
	 * The sum of matrix r has to be the number of examinees multiplied by number of items N x p.
	 * @param r the r matrix.
	 * @param N number of examinees.
	 * @param p number of items.
	 * @return true if NxP is equals to sum of matrix r
	 */
	bool test_r ( std::vector<matrix<double> > &r, int N, int p );

	/**
	 * Overloading test_r function.
	 * @param r the r matrix.
	 * @param N number of examinees.
	 * @param p number of items.
	 * @return true if NxP is equals to sum of matrix r
	 * @see test_r
	 */
	bool test_r ( matrix<double> &r, int N, int p );
}

#endif /* TEST_TEST_H_ */
