/*
 * gpc.h
 *
 *  Created on: Sep 26, 2016
 *      Author: Milder
 */

#ifndef POLYTOMOUS_MODEL_GPC_H_
#define POLYTOMOUS_MODEL_GPC_H_

#include "model.h"

namespace latentregpp {

namespace polytomous {

class gpc: public polytomous::model {
public:
	gpc();
	virtual ~gpc();
	gpc(int, std::vector<int>*);

	double exp_sum(std::vector<double>&, const optimizer_vector&, int k);

	double Pik(std::vector<double>&, const optimizer_vector&, int i, int k);
};

}

} /* namespace latentregpp */


#endif /* POLYTOMOUS_MODEL_GPC_H_ */
