#ifndef SIMANN_HPP
#define SIMANN_HPP

#include <Optimization/Optim.hpp>
#include <Optimization/Heuristic/Parameter.hpp>
#include <cmath>
#include <random>
#include <time.h>
#include <math.h>

namespace spgo
{
  template <typename Function, typename Info>
  class SimulatedAnnealing : public Optim<Function, Parameter, Info>
  {
    public:
      SimulatedAnnealing()
      {
        initialice(0.0001, 100);
      }

      SimulatedAnnealing(double tol, double rad)
      {
        initialice(tol, rad);
      }

      ~SimulatedAnnealing(){}
            
      int run(std::vector<Function> f, Parameter &p, Info info)
      {
        std::random_device         rd;
        std::mt19937               e2(rd());
        std::normal_distribution<> dist(-.5, 1);
        Parameter*                 new_p;
        double                     eval,
                                   new_eval,
                                   change;

        change = 1000;
        eval   = (f[0])(&p, info)[0];
        new_eval = eval + 100;
                
        while(change > tol)
        {
          new_p = p.getNeighbor(rad);
          new_eval = (f[0])(new_p, info)[0];

          if(new_eval < eval)
          {
            change = abs(eval - new_eval);
            delete p.parameters[0];
            delete p.parameters[1];
            p.parameters[0] = new_p->parameters[0];
            p.parameters[1] = new_p->parameters[1];
            eval = new_eval;
          }
          else
          {
            delete new_p->parameters[0];
            delete new_p->parameters[1];
          }

          delete new_p;

          rad = rad <= 1 ? rad + 1 : rad + dist(e2);
        }

        return 0;
      }

    private:
      double tol,
             rad;

      void initialice(double tol, double rad)
      {
        this->tol = tol;
        this->rad = rad;
      }
  };
}

#endif