/**
 * estep.cpp
 *
 *  Created on: 13/04/2016
 *      Author: Milder
 */

#include "../../dichotomous/estimation/estep.h"

namespace latentregpp {

namespace dichotomous {

void Estep ( estimation_data &data, int current ) {
	//Number of items
	int &p = data.p;
	//Number of response patterns
	int &s = data.s;
	//Number of quadrature points
	int &G = data.G;
	//Matrix of response patterns
	matrix<char> &Y = data.Y;
	//Frequency of each pattern
	std::vector<int> &nl = data.nl;
	//Latent trait vectors
	matrix<double> &theta = data.theta;
	//Weights
	std::vector<double> &w = data.w;
	//Vector of parameters of the items
	std::vector<optimizer_vector> &zeta = data.zeta[current];
	//f
	std::vector<double> &f = data.f;
	f.assign(f.size(), 0);

	//Matrix correct that has been answered correctly
	matrix<int> &correct = data.correct;

	//pi matrix
	matrix<double> &pi = data.pi;


	// Probability matrix P
	matrix<double> &P = data.P;

	//r matrix
	matrix<double> &r = data.r;
	r.reset();

	/**
	 * Computing each element of matrix P
	 * P_gi
	 */
	#pragma omp parallel for schedule(dynamic)
	for ( int g = 0; g < G; ++g ) {
		std::vector<double> &theta_g = *theta.get_pointer_row(g);
		for ( int i = 0; i < p; ++i ) {
			P(g, i) = data.m->P(theta_g, zeta[i]);
		}
	}

	double integral_l = 0;

	/**
	 * Calcule the Pi_gl. Equation (7) from IRT_engineers document
	 * this has a complexity O(G*s*p)
	 */

	//Patterns
	#pragma omp parallel for schedule(dynamic) reduction(+:integral_l)
	for ( int l = 0; l < s; ++l ) {
		integral_l = 0;
		//Quadrature points
		for ( int g = 0; g < G; ++g ) {
			double &pi_gl = pi(g, l) = w[g];
			//Items
			for ( int i = 0; i < p; ++i )
				pi_gl *= Y(l, i) ? P(g, i) : 1 - P(g, i);
			/**
			 * As denominator for a response pattern l is the summation over the latent traits
			 * here pi(g, l) is added to integral_l
			 * */
			integral_l += pi_gl;
		}

		for ( int g = 0; g < G; ++g ) {
			double &pi_gl = pi(g, l);
			pi_gl /= integral_l;
		}
	}

	#pragma omp parallel for schedule(dynamic)
	for ( int g = 0; g < G; ++g ) {
		for ( int l = 0; l < s; ++l ) {
			double &pi_gl = pi(g, l);
			f[g] += nl[l] * pi_gl;
			for ( int i = 0; i < correct.columns(l); ++i )
				r(g, correct(l, i)) += nl[l] * pi_gl;
		}
	}

//	//Asserting pi correctness
//	bool pi_ok = test_pi(pi);
//	assert(("Each column of pi matrix must sum 1.0", pi_ok));
}

}

} /* namespace latentregpp */
