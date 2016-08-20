#ifndef PARAMETER_HPP
#define PARAMETER_HPP

#include <Types/Generic.hpp>
#include <vector>

namespace spgo
{
  class Parameter
  {
    public:
      Parameter(){}

      virtual ~Parameter()
      {
        for(int i = 0; i < parameters.size(); i++)
        {
           delete parameters[i];
        }
      }

      virtual void copy(Parameter* p) = 0;

      virtual Parameter* getNeighbor(double rad) = 0;

      std::vector<Any*> parameters;
  };
}

#endif