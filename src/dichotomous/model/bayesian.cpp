/*
 * bayesian.cpp
 *
 *  Created on: 26/10/2016
 *      Author: jhonatan
 */

#include "bayesian.h"

namespace latentregpp {

namespace dichotomous {

bayesian::bayesian() : model(model_type::bayesian, THREE_PARAMETERS) {}

bayesian::bayesian(matrix<double> c) : model(model_type::bayesian, THREE_PARAMETERS) {
	this->c_values = c;
    //this->type = model_type::bayesian;
    //this->parameters = THREE_PARAMETERS;
    //Rprintf("%f",c_values(2,3));
    //Rprintf("I am in construct\n");
    /*for(int p = 0; p < c_values.rows(); p++)
    {
        for(int q = 0; q < c_values.columns(p); q++) {
            Rprintf("%f ",c_values(p,q));
        }
        Rprintf("\n");
    }*/
}

bayesian::bayesian(std::vector<optimizer_vector> c) {
    //this->c_values = c;
    Rprintf("Cosntructor with optimizer_vector");
}

double bayesian::P(std::vector<double> &theta, const optimizer_vector &parameters, int i) {
	double gamma_parameter = parameters(parameters.size() - 1);

	//uncommented line below for reparameter a c value [0,1] from gamma in R
	//double c = 1.0 / (1.0 + exp(-gamma_parameter));

	//need to put c initial value here
	int item = i;
	//double c = 0.1;
	double c = c_values(item,c_values.columns(item)-1);
	c = 1.0 / (1.0 + exp(-c));

	//Initialized with gamma value
	double eta = parameters(parameters.size() - 1);

	//Computing dot product
	for ( size_t j = 0; j < theta.size(); ++j )
		eta += parameters(j) * theta[j];

	/**three different formulas**/

	//with clear
	//double P = (1.0 / (1.0 + exp(-gamma_parameter))) + (1.0 - (1.0 / (1.0 + exp(-gamma_parameter)))) / (1.0 + exp(-eta));

	//without clear
	//double P = (1.0 / (1.0 + std::exp(-gamma_parameter))) + (1.0 / (1.0 + std::exp(gamma_parameter))) * (1.0 / (1.0 + std::exp(-eta)));

	//with reparameter
	double P = c + (1.0 - c) / (1.0 + exp(-eta));

	P = std::max(P, LOWER_BOUND_);
	P = std::min(P, UPPER_BOUND_);
	return P;
}

void bayesian::set_c_values(std::vector<optimizer_vector> c) {
    int row = c.size();
    int col = c[0].size();
    c_values = matrix<double>(row,col);
    for(int p = 0; p < c.size() ;++p) {
        for(int j = 0; j < c[p].size() ;++j) {
            c_values(p,j) = c[p](j);
        }
    }
    
    /*Rprintf("Print values for c_values matrix atributte");
    for(int x = 0; x < c_values.rows() ;++x) {
        for(int y = 0; y < c_values.columns(x) ;++y) {
            Rprintf("%f ",c_values(x,y));
        }
        Rprintf("\n");
    }*/
}

bayesian::~bayesian() {}

}

} /* namespace latentregpp */

