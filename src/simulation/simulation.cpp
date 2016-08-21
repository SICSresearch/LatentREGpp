/*
 * simulation.cpp
 *
 *  Created on: 1/06/2016
 *      Author: Milder
 */

#include "simulation.h"
#define START_CLOCK		double start = omp_get_wtime();
#define END_CLOCK		double stop = omp_get_wtime();
#define REPORT_TIME     double elapsed = stop - start;\
						std::cout << "Time elapsed: " << elapsed << "s." << '\n';


namespace irtpp {

simulation::simulation() {

}

void simulation::simulate ( int model, int d, int start, int end, std::string folder,
							std::string name, double dif, bool dicho,
							std::string quadrature_technique, int G, std::vector<int> cluster,
							std::string custom_initial_values_filename ) {

	std::ofstream report_parameters;
	std::stringstream ss;
	ss << folder << "/estimation-" << name << '-' << start << '-' << end;
	if ( quadrature_technique == QMCEM ) ss << "-G=" << G;
	ss << (dicho ? "-dicho-package" : "-poly-package") << ".csv";

	std::string parameters = ss.str();
	report_parameters.open(parameters.c_str());
	report_parameters.precision(4);

	ss.str("");
	ss << folder << "/" << name;
	const std::string base_name = ss.str();

	ss.str("");
	ss << folder << "/INI/INITIAL-" << name;
	const std::string initial_base_name = ss.str();

	for ( int i = start; i <= end; ++i ) {
		matrix<char> Y;
		input<char> in;

		ss.str("");
		ss << base_name << i << ".csv";

		std::string file_name = ss.str();
		in.import_data(file_name, Y);

		ss << " imported. Running with " << (dicho ? "dichotomous" : "polytomus") << " package";
		if ( quadrature_technique == QMCEM )
			ss << " and Sobol G=" << G;
		else
			ss << " and Gaussian";

		std::cout << ss.str() << std::endl;

		std::string initial_values = custom_initial_values_filename;
		if ( custom_initial_values_filename == BUILD ) {
			ss.str("");
			ss << initial_base_name << i << ".csv";
			initial_values = ss.str();
		}

		std::cout << "Initial values from file: " << initial_values << std::endl;

		if ( dicho ) {
			START_CLOCK
			dichotomous::estimation e(Y, d, model, dif, cluster, quadrature_technique, G, EMPTY_INTEGER_VECTOR, initial_values);
			e.EMAlgortihm();

			END_CLOCK
			REPORT_TIME
			e.print_item_parameters(report_parameters, elapsed);
		} else {
			START_CLOCK
			polytomous::estimation e(Y, d, model, dif, cluster, quadrature_technique, G, EMPTY_INTEGER_VECTOR, initial_values);
			e.EMAlgortihm();

			END_CLOCK
			REPORT_TIME
			e.print_item_parameters(report_parameters, elapsed);
		}
	}

	report_parameters.close();
}

void simulation::simulate ( int model, int d, int iterations, std::string folder,
							std::string name, int interval, double dif, bool dicho,
							std::string quadrature_technique, int G,
							std::vector<int> cluster, std::string custom_initial_values_filename) {
	for ( int i = 1; i <= iterations; i += interval ) {
		simulate(model, d, i, i + interval - 1, folder, name, dif, dicho,
				 quadrature_technique, G, cluster, custom_initial_values_filename);
	}
}

void simulation::run_single ( int model, int d, std::string filename, double dif, bool dicho,
							   std::string quadrature_technique, int G, std::vector<int> cluster,
							   std::string custom_initial_values_filename ) {
	if ( dicho ) run_single_dichotomous(model, d, filename, dif, quadrature_technique, G, cluster, custom_initial_values_filename );
	else		 run_single_polytomous(model, d, filename, dif, quadrature_technique, G, cluster, custom_initial_values_filename );
}

void simulation::run_single_polytomous ( int model, int d, std::string filename, double dif,
										std::string quadrature_technique, int G, std::vector<int> cluster,
									    std::string custom_initial_values_filename ) {
	matrix<char> Y;
	input<char> in;
	in.import_data(filename, Y);
	std::cout << "Data imported from " << filename << std::endl;

	START_CLOCK

	polytomous::estimation e(Y, d, model, dif, cluster, quadrature_technique, G, EMPTY_INTEGER_VECTOR, custom_initial_values_filename);
	e.EMAlgortihm();

	END_CLOCK
	e.print_item_parameters();
	REPORT_TIME

	std::stringstream ss;
	ss << filename << "-estimation.csv";
	std::ofstream out(ss.str());
	e.print_item_parameters(out, elapsed);
	out.close();
}

void simulation::run_single_dichotomous ( int model, int d, std::string filename, double dif,
										  std::string quadrature_technique, int G, std::vector<int> cluster,
										  std::string custom_initial_values_filename ) {
	matrix<char> Y;
	input<char> in;
	in.import_data(filename, Y);
	std::cout << "Data imported from " << filename << std::endl;

	START_CLOCK

	dichotomous::estimation e(Y, d, model, dif, cluster, quadrature_technique, G, EMPTY_INTEGER_VECTOR, custom_initial_values_filename);
	e.EMAlgortihm();

	END_CLOCK
	e.print_item_parameters();
	REPORT_TIME

	std::stringstream ss;
	ss << filename << "-estimation.csv";
	std::ofstream out(ss.str());
	e.print_item_parameters(out, elapsed);
	out.close();
}

simulation::~simulation() {

}

} /* namespace irtpp */
