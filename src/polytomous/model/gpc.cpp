/*
 * gpc.cpp
 *
 *  Created on: Sep 26, 2016
 *      Author: Milder
 */

#include "gpc.h"

namespace latentregpp {

namespace polytomous {

gpc::gpc(int d, std::vector<int> *categories_item) :
	model(model_type::gpc, TWO_PARAMETERS, d, categories_item) {
	printf("GPC\n");
}

gpc::~gpc() {}

double gpc::exp_sum(std::vector<double>& theta, const optimizer_vector &parameters, int k) {
	int d = theta.size();
	double sum = 0;

	for (int h = 1; h <= k; ++h) {
		double eta = -parameters(d + h);
		//Dot product
		for ( int i = 0; i < d; ++i )
			eta += parameters(i) * theta[i];
		sum += eta;
	}

	return std::exp(sum);
}

double gpc::Pik(std::vector<double> &theta, const optimizer_vector &parameters, int i, int k) {
	int mi = (*categories_item)[i];

	printf("GPC\n");

	double sum = 0;
	for (int c = 1; c < mi; ++c) {
		sum += exp_sum(theta, parameters, c);
	}

	double num = ( k == 0 ? 1.0 : exp_sum(theta, parameters, k) );
	double den = (1.0 + sum);

	printf("%d %d: %.5lf %.5lf\n", i + 1, k + 1, num, den);

	double P_ik = num / den;

	P_ik = std::max(P_ik, LOWER_BOUND_);
	P_ik = std::min(P_ik, UPPER_BOUND_);

	return P_ik;
}

}

} /* namespace latentregpp */
