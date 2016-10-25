/*
 * estimation.cpp
 *
 *  Created on: 13/04/2016
 *      Author: Milder
 */

#include "../../dichotomous/estimation/estimation.h"

namespace latentregpp {

namespace dichotomous {

estimation::estimation(matrix<char> &dataset, unsigned int d, int themodel,
					   double convergence_difference,
					   matrix<double> theta,
					   std::vector<double> weights,
					   std::vector<int> individual_weights,
					   std::vector<int> clusters,
					   matrix<double> initial_values ) {
	/**
	 * Object to allocate all data needed in estimation process
	 * */
	data = estimation_data(d);
	data.dataset = &dataset;

	//-------------------------------------------------------------------------------------

	//Number of examinees
	int &N = data.N;

	//Number of items
	int &p = data.p;

	//Number of response patterns (s <= N)
	int &s = data.s;

	//Matrix of response patterns. Its size is s x p
	matrix<char> &Y = data.Y;

	//Frequency of each pattern
	std::vector<int> &nl = data.nl;

	//Matrix correct that has been answered correctly
	matrix<int> &correct = data.correct;

	//Matrix of response patterns and their frequency
	std::map<std::vector<char>, std::vector<int> > &patterns = data.patterns;

	double &loglikelihood = data.loglikelihood;

	//-------------------------------------------------------------------------------------


	for ( int i = 0; i < dataset.rows(); ++i )
		patterns[dataset.get_row(i)].push_back(i);

	Y = matrix<char>();
	nl = std::vector<int>(patterns.size());

	if ( individual_weights.empty() ) individual_weights = std::vector<int>(patterns.size(), 1);

	int l = 0;
	for ( auto it : patterns ) {
		Y.add_row(it.first);
		nl[l] = it.second.size() * individual_weights[l];
		++l;
	}

	N = dataset.rows();
	s = Y.rows();
	p = Y.columns(0);

	correct = matrix<int>(s);
	for ( int l = 0; l < s; ++l )
		for ( int i = 0; i < p; ++i )
			if ( Y(l, i) )
				correct.add_element(l, i);

	data.theta = theta;
	data.w = weights;
	data.G = theta.rows();
	build_matrixes();

	switch ( themodel ) {
		case model_type::onepl:
			data.m = new onepl();
			break;
		case model_type::twopl:
			data.m = new twopl();
			break;
		case model_type::threepl:
			data.m = new threepl();
			break;
		default:
			data.m = new twopl();
	}


	if ( d == 1 ) compute_1D_initial_values();
	else {
		//Pinned items in multidimensional case (the first of each dimension)
		std::set<int> &pinned_items = data.pinned_items;

		//clusters size MUST be equal to the number of dimensions
		if ( clusters.size() == d ) {
			int pinned = 0;
			for ( unsigned int i = 0; i < clusters.size(); ++i ) {
				pinned_items.insert(pinned);
				pinned += clusters[i];
			}
		}

		if ( initial_values.rows() > 0 )
			load_multi_initial_values(initial_values);
	}

	//Configurations for the estimation
	loglikelihood = NOT_COMPUTED;
	this->convergence_difference = convergence_difference;
	this->iterations = 0;	
}

void estimation::build_matrixes() {
	//Number of items
	int &p = data.p;
	//Number of response patterns (s <= N)
	int &s = data.s;
	//Number of quadrature points
	int &G = data.G;
	//Matrix r. Needed in Estep and Mstep
	matrix<double> &r = data.r;

	/**
	 * Probability matrix P
	 *
	 * P_gi
	 *
	 * P_gi means the probability that an individual has selected the correct answer
	 *
	 *
	 * The purpose of this matrix is to allocate the value of P_gi
	 * to avoid recompute them while matrix Pi and r are computed in EStep
	 * */
	matrix<double> &P = data.P;

	//Matrix pi
	matrix<double> &pi = data.pi;
	//f
	std::vector<double> &f = data.f;

	//Builds r and P matrixes
	P = matrix<double>(G, p);
	pi = matrix<double>(G, s);
	r = matrix<double>(G, p);
	f = std::vector<double>(G);
}


void estimation::load_multi_initial_values ( matrix<double> &mt ) {
	//Dimension
	int &d = data.d;
	//Parameters of the items
	std::vector<optimizer_vector> &zeta = data.zeta[0];
	//Number of items
	int &p = data.p;

	zeta = std::vector<optimizer_vector>(p);
	int total_parameters = data.m->parameters == ONE_PARAMETER ? 1 : data.m->parameters - 1 + d;
	
	for ( int i = 0; i < p; ++i ) {
		zeta[i] = optimizer_vector(total_parameters);
		if ( data.m->parameters == ONE_PARAMETER ) {
			zeta[i](0) = mt(i, d);
		} else {
			for ( int j = 0; j < total_parameters; ++j )
				zeta[i](j) = mt(i, j);
			if ( data.m->parameters == THREE_PARAMETERS ) {
				double &c = zeta[i](total_parameters - 1);
				c = std::log(c / (1.0 - c));
			}
		}
	}

	//Items that will not be estimated
	std::set<int> &pinned_items = data.pinned_items;

	if ( d > 1 && pinned_items.empty() ) {
		int items_for_dimension = (p + d - 1) / d;
		for ( int i = 0, j = 0; i < p; i += items_for_dimension, ++j )
			pinned_items.insert(i);
	}

	int j = 0;
	for ( auto pinned : pinned_items ) {
		optimizer_vector &item = zeta[pinned];
		for ( int h = 0; h < d; ++h )
			if ( h != j )
				item(h) = 0;
		++j;
	}

	data.loglikelihood = NOT_COMPUTED;
	iterations = 0;
}

void estimation::compute_1D_initial_values() {
	//Parameters of the items
	std::vector<optimizer_vector> &zeta = data.zeta[0];
	//Number of examinees
	int &p = data.p;
	//Matrix of answers of the examinees
	matrix<char> &dataset = *data.dataset;

	zeta = std::vector<optimizer_vector>(p);
	int total_parameters = data.m->parameters;

	for ( int i = 0; i < p; ++i ) {
		zeta[i] = optimizer_vector(total_parameters);
		for ( int j = 0; j < total_parameters; ++j )
			zeta[i](j) = DEFAULT_INITIAL_VALUE;
	}

	std::vector<double> alpha, gamma;
	find_initial_values(dataset, alpha, gamma);

	for ( int i = 0; i < p; ++i ) {
		optimizer_vector &item_i = zeta[i];

		if ( data.m->parameters > ONE_PARAMETER ) {
			item_i(0) = alpha[i];
			item_i(1) = gamma[i];
			if ( data.m->parameters == THREE_PARAMETERS ) item_i(2) = DEFAULT_C_INITIAL_VALUE;
		} else {
			item_i(0) = gamma[i];
		}
	}

	data.loglikelihood = NOT_COMPUTED;
}

void estimation::EMAlgorithm ( bool verbose ) {
	if ( verbose ) Rprintf("EM Algorithm started\n");
	double dif = 0.0;
	iterations = 0;
	int current;
	do {
		current = iterations % ACCELERATION_PERIOD;
		if ( current == 2 )
			ramsay(data.zeta, data.pinned_items);

		Estep(data, current);
		dif = Mstep(data, current);
		++iterations;
		if ( verbose ) Rprintf("\rIteration: %u \tMax-Change: %.6lf", iterations, dif);
	} while ( dif >= convergence_difference && iterations < MAX_ITERATIONS );
	if ( verbose ) Rprintf("\n");
}

double estimation::log_likelihood() {
	double &loglikelihood = data.loglikelihood;
	if ( loglikelihood != NOT_COMPUTED ) return loglikelihood;

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
	std::vector<optimizer_vector> &zeta = data.zeta[iterations % ACCELERATION_PERIOD];

	// Probability matrix P
	matrix<double> &P = data.P;

	#pragma omp parallel for schedule(dynamic)
	for ( int g = 0; g < G; ++g ) {
		std::vector<double> &theta_g = *theta.get_pointer_row(g);
		for ( int i = 0; i < p; ++i ) {
			P(g, i) = data.m->P(theta_g, zeta[i]);
		}
	}

	double integral_l = 0;
	double loglikelihood_temp = 0;

	//Patterns
	#pragma omp parallel for schedule(dynamic) reduction(+:integral_l,loglikelihood_temp)
	for ( int l = 0; l < s; ++l ) {
		integral_l = 0;
		//Quadrature points
		for ( int g = 0; g < G; ++g ) {
			double pi_gl = w[g];
			//Items
			for ( int i = 0; i < p; ++i )
				pi_gl *= Y(l, i) ? P(g, i) : 1 - P(g, i);
			integral_l += pi_gl;
		}
		loglikelihood_temp += nl[l] * std::log(integral_l);
	}

	return loglikelihood = loglikelihood_temp;
}

void estimation::EAP ( bool all_factors ) {
	//Number of response patterns
	int &s = data.s;
	//Number of quadrature points
	int &G = data.G;
	//Dimension
	int &d = data.d;
	//Latent trait vectors
	matrix<double> &theta = data.theta;
	//pi matrix
	matrix<double> &pi = data.pi;

	Estep(data, iterations % ACCELERATION_PERIOD);

	//Patterns
	std::map<std::vector<char>, std::vector<int> > &patterns = data.patterns;
	//Latent traits
	std::vector<optimizer_vector> &latent_traits = data.latent_traits;
	latent_traits = std::vector<optimizer_vector>(s);

	int l = 0;
	for ( auto pt : patterns ) {
		optimizer_vector &theta_l = latent_traits[l];
		theta_l = optimizer_vector(d);

		for ( int g = 0; g < G; ++g ) {
			std::vector<double> theta_g = *theta.get_pointer_row(g);
			for ( int h = 0; h < d; ++h )
				theta_l(h) += theta_g[h] * pi(g, l);
		}
		++l;
	}

	if ( all_factors ) latent_traits_by_individuals();
}


void estimation::MAP ( bool all_factors ) {
	std::vector<optimizer_vector> &latent_traits = data.latent_traits;
	int &s = data.s;

	int current_zeta = iterations % ACCELERATION_PERIOD;
	for ( int l = 0; l < s; ++l ) {
		dlib::find_max_using_approximate_derivatives(dlib::bfgs_search_strategy(),
							   dlib::objective_delta_stop_strategy(OPTIMIZER_DELTA_STOP),
							   posterior(l, current_zeta, &data), latent_traits[l], -1);
	}	

	if ( all_factors ) latent_traits_by_individuals();
}



estimation::posterior::posterior (int l, int cur, estimation_data *d) :
		l(l), current_zeta(cur), data(d) { }

double estimation::posterior::operator() ( const optimizer_vector& theta_l ) const {
	int p = data->p;
	int d = data->d;
	matrix<char> &Y = data->Y;
	std::vector<optimizer_vector> &zeta = data->zeta[current_zeta];

	double value = 0.0;
	for ( int h = 0; h < d; ++h )
		value += theta_l(h) * theta_l(h);

	value = std::exp(-0.5 * value);
	value *= p;

	for ( int i = 0; i < p; ++i )
		value += Y(l, i) ? std::log(data->m->P(theta_l, zeta[i])) : std::log(1 - data->m->P(theta_l, zeta[i]));

	return value;
}

void estimation::latent_traits_by_individuals () {
	std::vector<optimizer_vector> &latent_traits = data.latent_traits;
	std::vector<optimizer_vector> temp_latent_traits = latent_traits;
	latent_traits = std::vector<optimizer_vector>(data.N);
	std::map<std::vector<char>, std::vector<int> > &patterns = data.patterns;
	int l = 0;
	for ( auto pt : patterns ) {
		for ( size_t j = 0; j < pt.second.size(); ++j )
			latent_traits[pt.second[j]] = temp_latent_traits[l];
		++l;
	}
}

unsigned int estimation::get_iterations ( ) {
	return iterations;
}

estimation::~estimation() {

}

} /* namespace dichotomous */

} /* namespace latentregpp */
