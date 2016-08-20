#ifndef OPTIM_HPP
#define OPTIM_HPP
#include <vector>

namespace spgo
{

enum STATE{CONVERGENCE, NO_CONVEGENCE, INVALID};

template <typename Function, typename Parameter, typename Info>
class Optim
{
	public:
		
		Optim(){}
		
		virtual int run(std::vector<Function> f, Parameter &p, Info info) = 0;
		
		virtual ~Optim(){};
};

}

#endif