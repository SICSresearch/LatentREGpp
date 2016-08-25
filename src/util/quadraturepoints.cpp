/**
 * quadraturepoints.cpp
 *
 *  Created on: 22/04/2016
 *      Author: Milder
 */

#include "quadraturepoints.h"
#include "input.h"
#include "matrix.h"

namespace latentregpp {

	void compute_quadrature_points ( std::vector<double> &q, int r,
								std::vector<double> &current_trait,
				  	  	  	    matrix<double> &latent_trait ) {
		if ( r == 0 ) {
			latent_trait.add_row(current_trait);
			return;
		}
		for ( size_t i = 0; i < q.size(); ++i ) {
			current_trait.push_back(q[i]);
			compute_quadrature_points(q, r - 1, current_trait, latent_trait);
			current_trait.pop_back();
		}
	}

	void compute_and_save_quadrature_points(int G, int d) {
		std::vector<double> q;
		input<double> in;
		//Quadrature points loaded from file
		std::stringstream ss;
		ss << G;
		std::string filename = getwd() + "data/quadrature" + ss.str() + "_in.data";
		if ( !in.import_data(filename, q) ) {
			std::cout << "Your filename " << filename << '\n';
			std::cout << "Filename must have this form: " << '\n';
			std::cout << "quadratureG_in.data and be allocated in /data/ \n";
			return;
		}

		double f = sqrt(2);
		for ( size_t i = 0; i < q.size(); ++i )
			q[i] *= f;


		/**
		 * Here the latent trait vectors are computed
		 * using a backtracking approach
		 * */

		//latent trait vectors
		matrix<double> latent_trait;
		//current latent trait vector, it starts void
		std::vector<double> current_trait;
		//computing latent trait vector for the given dimension
		compute_quadrature_points(q, d, current_trait, latent_trait);



		std::cout << "Latent trait vectors. Total[" << latent_trait.rows() << "]\n";
		for ( int i = 0; i < latent_trait.rows(); ++i ) {
			std::cout << '[';
			for ( int j = 0; j < latent_trait.columns(i); ++j ) {
				if (j) std::cout << ", ";
				std::cout << latent_trait(i, j);
			}
			std::cout << "]\n";
		}

		std::ofstream out;
		ss.str("");
		ss << d;
		filename = getwd() + "data/quadrature" + ss.str() + "_computed.data";
		out.open(filename.c_str());
		for ( int i = 0; i < latent_trait.rows(); ++i ) {
			for ( int j = 0; j < latent_trait.columns(i); ++j )
				out << latent_trait(i, j) << ' ';
			out << '\n';
		}
		out.close();
	}

	void compute_weights ( std::vector<double> &w, int r,
						   double current_weight, std::vector<double> &weights ) {
		if ( r == 0 ) {
			weights.push_back(current_weight);
			return;
		}
		for ( size_t i = 0; i < w.size(); ++i ) {
			current_weight *= w[i];
			compute_weights(w, r - 1, current_weight, weights);
			current_weight /= w[i];
		}
	}

	void compute_and_save_weights ( int G, int d ) {
		std::vector<double> w;
		input<double> in;
		std::stringstream ss;
		ss << G;
		std::string filename = getwd() + "data/weights" + ss.str() + "_in.data";
		if ( !in.import_data(filename, w) ) {
			std::cout << "Your filename " << filename << '\n';
			std::cout << "Filename must have this form: " << '\n';
			std::cout << "weightsG_in.data and be allocated in /data/\n";
			return;
		}

		double pi = acos(-1);
		double f = pow(sqrt(pi), d);
		for ( size_t i = 0; i < w.size(); ++i )
			w[i] /= f;

		/**
		 * Weights are computed as the same way than latent trait vectors
		 * Using backtracking
		 * */

		std::vector<double> weights;
		double current_weight = 1;
		compute_weights(w, d, current_weight, weights);

		std::cout << "Weights. Total[" << weights.size() << "]\n";
		for ( size_t i = 0; i < weights.size(); ++i ) {
			std::cout << weights[i] << '\n';
		}

		std::ofstream out;
		ss.str("");
		ss << d;
		filename = getwd() + "data/weights" + ss.str() + "_computed.data";
		out.open(filename.c_str());
		for ( size_t i = 0; i < weights.size(); ++i )
			out << weights[i] << '\n';
		out.close();
	}

	matrix<double> load_quadrature_points ( int d ) {
		std::stringstream ss;
		ss << d;
		std::string filename = getwd() + "data/quadrature" + ss.str() + "_computed.data";

		matrix<double> m;
		input<double> in;
		in.import_data(filename, m);
		return m;
	}

	std::vector<double> load_weights ( int d ) {
		std::stringstream ss;
		ss << d;
		std::string filename = getwd() + "data/weights" + ss.str() + "_computed.data";

		std::vector<double> v;
		input<double> in;
		in.import_data(filename, v);
		return v;
	}
}


