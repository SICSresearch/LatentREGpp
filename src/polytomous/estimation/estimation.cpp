/*
 * estimation.cpp
 *
 *  Created on: 13/04/2016
 *      Author: Milder
 */

#include "estimation.h"

namespace irtpp {

namespace polytomous {

estimation::estimation(matrix<char> &dataset, unsigned int d, int themodel,
					   double convergence_difference,
					   std::vector<int> number_of_items,
					   std::string quadrature_technique,
					   int quadrature_points,
					   std::vector<int> individual_weights,
					   std::string custom_initial_values_filename ) {
	/**
	 * Object to allocate all data needed in estimation process
	 * */
	data = estimation_data(d);
	data.dataset = &dataset;

	//-------------------------------------------------------------------------------------

	//Model to be used
	model &m = data.m;

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

	//Number of categories by item
	std::vector<int> &categories_item = data.categories_item;

	//Matrix of response patterns and their frequency
	std::map<std::vector<char>, std::vector<int> > &patterns = data.patterns;

	double &loglikelihood = data.loglikelihood;

	//-------------------------------------------------------------------------------------

	bool dichotomous = false;
	for ( int i = 0; i < dataset.rows(); ++i ) {
		for ( int j = 0; j < dataset.columns(i); ++j )
			dichotomous |= dataset(i, j) == 0;
	}

	//If dataset is dichotomous, it's converted to polytomous with two categories
	if ( dichotomous ) {
		for ( int i = 0; i < dataset.rows(); ++i )
			for ( int j = 0; j < dataset.columns(i); ++j )
				++dataset(i, j);
	}


	//Matrix of response patterns and their frequency
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

	//Number of categories of each item
	categories_item = std::vector<int>(p);
	for ( int i = 0; i < p; ++i ) {
		int max_category = -1;
		for ( int l = 0; l < s; ++l ) {
			//Number of categories of an item is defined as the max category found in the answers
			if ( Y(l, i) > max_category )
				max_category = Y(l, i);
			//Checking if data is dichotomous
		}
		categories_item[i] = max_category;
	}

	if ( quadrature_technique == QMCEM )
		sobol_quadrature(quadrature_points);
	else
		gaussian_quadrature();

	//Pinned items in multidimensional case (the first of each dimension)
	std::set<int> &pinned_items = data.pinned_items;

	//Number of items size MUST be equal to the number of dimensions
	if ( number_of_items.size() == d ) {
		int pinned = 0;
		for ( unsigned int i = 0; i < number_of_items.size(); ++i ) {
			pinned_items.insert(pinned);
			pinned += number_of_items[i];
		}
	}

	//Configurations for the estimation
	loglikelihood = NOT_COMPUTED;
	m = model(themodel, d, &categories_item);
	this->convergence_difference = convergence_difference;
	this->iterations = 0;
	this->custom_initial_values_filename = custom_initial_values_filename;
}



void estimation::build_matrixes() {
	//Number of items
	int &p = data.p;

	//Number of response patterns (s <= N)
	int &s = data.s;

	//Number of quadrature points
	int &G = data.G;

	//Number of categories by item
	std::vector<int> &categories_item = data.categories_item;

	//Matrix r. Needed in Estep and Mstep
	std::vector<matrix<double> > &r = data.r;

	/**
	 * Probability matrix P
	 *
	 * P_gik
	 *
	 * P_gik means the probability that an individual has selected the category k
	 * to item i and belongs to group g
	 *
	 *
	 * The purpose of this matrix is to allocate the value of P_gik
	 * to avoid recompute them while numerators and denominators in Estep are computed
	 * */
	std::vector<matrix<double> > &P = data.P;

	/**
	 * Matrix of probabilities pi, denominators vector and matrix of numerators
	 * needed in Estep
	 * */
	matrix<double> &pi = data.pi;

	//Builds r and P matrixes
	r = std::vector<matrix<double> >(G);
	P = std::vector<matrix<double> >(G);
	for ( int g = 0; g < G; ++g ) {
		r[g] = matrix<double>();
		P[g] = matrix<double>();
		for ( int i = 0; i < p; ++i ) {
			r[g].add_row(categories_item[i]);
			P[g].add_row(categories_item[i]);
		}
	}
	pi = matrix<double>(G, s);
}


void estimation::sobol_quadrature (int g) {
	//Dimension
	int &d = data.d;

	//Number of quadrature points
	int &G = data.G;

	//Latent trait vectors
	matrix<double> &theta = data.theta;

	//Weights
	std::vector<double> &w = data.w;

	input<double> in;
	std::stringstream ss;
	ss << "data/sobol" << d << ".data";
	in.import_data(ss.str(), theta);

	G = g;

	w = std::vector<double>(G, 1.0);
	build_matrixes();
}

void estimation::gaussian_quadrature () {
	//Dimension
	int &d = data.d;

	//Number of quadrature points
	int &G = data.G;

	//Latent trait vectors
	matrix<double> &theta = data.theta;

	//Weights
	std::vector<double> &w = data.w;

	/**
	 * Number of quadrature points (G) is computed based on
	 * MAX_NUMBER_OF_QUADRATURE_POINTS and dimension of the problem, in this way
	 *
	 *
	 * G will be in 1dimension = 40 ---> 40^1 = 40
	 * 				2dimension = 20 ---> 20^2 = 400
	 * 				3dimension = 10 ---> 10^3 = 1000
	 * 				> 4dimension = 5 ---> 5^d
	 * */

	// Latent trait vectors loaded from file
	theta = load_quadrature_points(d);

	// Weights loaded from file
	w = load_weights(d);

	G = theta.rows();
	build_matrixes();
}

void estimation::load_initial_values ( std::string filename ) {
	matrix<double> mt;
	input<double> in;
	in.import_data(filename, mt);

	//Dimension
	int &d = data.d;
	//Parameters of the items
	std::vector<optimizer_vector> &zeta = data.zeta[0];
	//Number of items
	int &p = data.p;
	//Model used in the problem
	model &m = data.m;
	//Number of categories of each item
	std::vector<int> &categories_item = data.categories_item;

	zeta = std::vector<optimizer_vector>(p);

	for ( int i = 0; i < p; ++i ) {
		int total_parameters = m.parameters == ONEPL ? categories_item[i] - 1 : categories_item[i] - 1 + d;
		zeta[i] = optimizer_vector(total_parameters);
		for ( int j = 0; j < total_parameters; ++j )
			zeta[i](j) = mt(i, j);
	}

	//Items that will not be estimated
	std::set<int> &pinned_items = data.pinned_items;

	/**
	 * It is supposed that there are p / d items for each dimension
	 * if the user does not specify them
	 *
	 *
	 * */

	if ( pinned_items.empty() ) {
		int items_for_dimension = (p + d - 1) / d;
		for ( int i = 0, j = 0; i < p; i += items_for_dimension, ++j )
			pinned_items.insert(i);
	}

	data.loglikelihood = NOT_COMPUTED;
	iterations = 0;
}


void estimation::initial_values() {
	//Parameters of the items
	std::vector<optimizer_vector> &zeta = data.zeta[0];
	//Dimension
	int &d = data.d;
	//Number of examinees
	int &N = data.N;
	//Number of items
	int &p = data.p;
	//Number of categories of each item
	std::vector<int> &categories_item = data.categories_item;
	//Model used in the problem
	model &m = data.m;
	//Matrix of answers of the examinees
	matrix<char> &dataset = *data.dataset;

	zeta = std::vector<optimizer_vector>(p);

	for ( int i = 0; i < p; ++i ) {
		int total_parameters = m.parameters == ONEPL ? categories_item[i] - 1 : categories_item[i] - 1 + d;
		zeta[i] = optimizer_vector(total_parameters);
		for ( int j = 0; j < total_parameters; ++j )
			zeta[i](j) = DEFAULT_INITIAL_VALUE;
	}

	if ( d == 1 ) {
		/**
		 * Here, it is necessary find dichotomous items for each polytomous item
		 * */

		for ( int i = 0; i < p; ++i ) {
			optimizer_vector &item_i = zeta[i];
			int mi = categories_item[i];

			matrix<char> data_dicho(N, mi - 1);
			for ( int k = 1; k < mi; ++k ) {
				for ( int j = 0; j < N; ++j )
					data_dicho(j, k - 1) = dataset(j, i) >= k + 1;
			}

			//std::cout << data << std::endl;
			//std::cout << data_dicho << std::endl;

			std::vector<double> alpha, gamma;
			find_initial_values(data_dicho, alpha, gamma);

			/**
			 * Real alpha for this item will be the average among all alphas computed
			 * */

			double a = mean(alpha);

			if ( m.parameters > ONEPL ) {
				item_i(0) = a;
				for ( int k = 0; k < mi - 1; ++k )
					item_i(k + 1) = gamma[k];
			} else {
				for ( int k = 0; k < mi - 1; ++k )
					item_i(k) = gamma[k];
			}
		}

		//Here one item (with maximum number of categories) is pinned
		int item_to_pin = 0, max_categories = 0;
		for ( int i = 0; i < p; ++i ) {
			if ( data.categories_item[i] > max_categories ) {
				max_categories = data.categories_item[i];
				item_to_pin = i;
			}
		}

		data.pinned_items.insert(item_to_pin);
	}

	else {

		//TODO Compute Alphas

		/**
		 * Multidimensional case
		 * */

		int alphas = m.parameters == TWOPL ? d : 0;
		/**
		 * Polytomous case
		 *
		 * Here, it is necessary find dichotomous items for each polytomous item
		 * */

		for ( int i = 0; i < p; ++i ) {
			optimizer_vector &item_i = zeta[i];
			int mi = categories_item[i];

			matrix<char> data_dicho(N, mi - 1);
			for ( int k = 1; k < mi; ++k ) {
				for ( int j = 0; j < N; ++j )
					data_dicho(j, k - 1) = dataset(j, i) >= k + 1;
			}

			std::vector<double> alpha, gamma;
			find_initial_values(data_dicho, alpha, gamma);

			//As there is more than one gamma, it is necessary iterate over the number of categories
			for ( int k = 0; k < mi - 1; ++k )
				item_i(alphas + k) = gamma[k];
		}

		//Items that will not be estimated
		std::set<int> &pinned_items = data.pinned_items;

		/**
		 * It is supposed that there are p / d items for each dimension
		 * if the user does not specify them
		 *
		 *
		 * */

		if ( pinned_items.empty() ) {
			int items_for_dimension = (p + d - 1) / d;
			for ( int i = 0, j = 0; i < p; i += items_for_dimension, ++j )
				pinned_items.insert(i);
		}

		int j = 0;
		for ( auto pinned : pinned_items ) {
			optimizer_vector &item = zeta[pinned];
			for ( int h = 0; h < d; ++h )
				item(h) = 0;
			item(j) = DEFAULT_INITIAL_VALUE;
			++j;
		}
	}

	data.loglikelihood = NOT_COMPUTED;
}
void estimation::EMAlgorithm() {
	if ( custom_initial_values_filename == NONE || custom_initial_values_filename == BUILD ) initial_values();
	else load_initial_values(custom_initial_values_filename);
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

		std::cout << "Iteration: " << iterations << " \tMax-Change: " << dif << std::endl;
	} while ( dif >= convergence_difference && iterations < MAX_ITERATIONS );
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
	//Model used
	model &m = data.m;
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
	//Number of categories
	std::vector<int> &categories_item = data.categories_item;

	// Probability matrix P
	std::vector<matrix<double> > &P = data.P;

	#pragma omp parallel for schedule(dynamic)
	for ( int g = 0; g < G; ++g ) {
		std::vector<double> &theta_g = *theta.get_pointer_row(g);
		for ( int i = 0; i < p; ++i ) {
			int mi = categories_item[i];
			for ( int k = 0; k < mi; ++k )
				P[g](i, k) = m.Pik(theta_g, zeta[i], i, k);
		}
	}

	double integral_l = 0;
	double loglikelihood_temp = 0;

	//Patterns
	#pragma omp parallel for schedule(dynamic) reduction(+:integral_l,loglikelihood_temp)
	for ( int l = 0; l < s; ++l ) {
		integral_l = 0;
		//Quadrature points
		for (int g = 0; g < G; ++g ) {
			double pi_gl = w[g];
			//Items
			for (int i = 0; i < p; ++i )
				pi_gl *= P[g](i, Y(l, i) - 1);
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

//	sobol_quadrature(10000);
//	build_matrixes();

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
	EAP(false);
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
	model &m = data->m;
	std::vector<optimizer_vector> &zeta = data->zeta[current_zeta];

	double value = 0.0;
	for ( int h = 0; h < d; ++h )
		value += theta_l(h) * theta_l(h);
	value = std::exp(-0.5 * value) / std::pow( std::sqrt(2.0 * PI_), d );

	for ( int i = 0; i < p; ++i )
		value += std::log(m.Pik(theta_l, zeta[i], i, Y(l, i) - 1));

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


void estimation::print_item_parameters ( ) {
	std::vector<optimizer_vector> &zeta = data.zeta[iterations % ACCELERATION_PERIOD];
	int &p = data.p;

	std::cout << "Finished after " << iterations << " iterations.\n";
	for ( int i = 0; i < p; ++i ) {
		std::cout << "Item " << i + 1 << '\n';
		for ( int j = 0; j < zeta[i].size(); ++j )
			std::cout << zeta[i](j) << ' ';
		std::cout << '\n';
	}
}

void estimation::print_item_parameters ( std::string filename, double elapsed ) {
	std::ofstream fout(filename);
	print_item_parameters(fout, elapsed);
	fout.close();
}

void estimation::print_item_parameters ( std::ofstream &fout, double elapsed ) {
	std::vector<optimizer_vector> &zeta = data.zeta[iterations % ACCELERATION_PERIOD];
	int &p = data.p;

	for ( int i = 0; i < p; ++i ) {
		for ( int j = 0; j < zeta[i].size(); ++j ) {
			if ( j ) fout << ';';
			fout << zeta[i](j);
		}
		fout << ';' << iterations << ';' << elapsed << '\n';
	}
}

void estimation::print_latent_traits ( ) {
	std::vector<optimizer_vector> &latent_traits = data.latent_traits;
	for ( size_t i = 0; i < latent_traits.size(); ++i ) {
		std::cout << i + 1 << "\t";
		for ( int j = 0; j < latent_traits[i].size(); ++j )
			std::cout << latent_traits[i](j) << ' ';
		std::cout << std::endl;
	}
}

void estimation::print_latent_traits ( std::string filename ) {
	std::ofstream out(filename);
	print_latent_traits(out);
	out.close();
}

void estimation::print_latent_traits ( std::ofstream &out ) {
	std::vector<optimizer_vector> &latent_traits = data.latent_traits;
	for ( size_t i = 0; i < latent_traits.size(); ++i ) {
		for ( int j = 0; j < latent_traits[i].size(); ++j ) {
			if ( j ) out << ';';
			out << latent_traits[i](j);
		}
		out << '\n';
	}
}

short estimation::get_iterations ( ) {
	return iterations;
}

estimation::~estimation() {}

}

} /* namespace irtpp */
