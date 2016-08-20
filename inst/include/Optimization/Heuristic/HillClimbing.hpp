#ifndef HC_HPP
#define HC_HPP

#include <Optimization/Optim.hpp>
#include <Optimization/Heuristic/Parameter.hpp>
#include <cmath>
#include <iostream>

namespace spgo
{
  template <typename Function, typename Info>
  class HillClimbing : public Optim<Function, Parameter, Info>
  {
    public:
      HillClimbing()
      {
        initialize(0.0001, 1);
      }

      HillClimbing(double tol, double rad)
      {
        initialize(tol, rad);
      }

      ~HillClimbing(){}

      int run(std::vector<Function> f, Parameter &p, Info info)
      {
        Parameter* new_p;
        double     eval;
        double     new_eval;
        int        iter;

        iter     = 0;
        eval     = (f[0])(&p, info)[0];
        new_eval = eval + 100;
                
        while(iter < 1000000)
        {
          new_p = p.getNeighbor(3);
          new_eval = (f[0])(new_p, info)[0];

          if(new_eval < eval)
          {
            p.copy(new_p);
            eval = new_eval;
          }

          delete new_p;
                    
          iter++;
        }

        return 0;
      }

      private:
        double tol;
        double rad;

        void initialize(double tol, double rad)
        {
          this->tol = tol;
          this->rad = rad;
        }
  };
}

#endif