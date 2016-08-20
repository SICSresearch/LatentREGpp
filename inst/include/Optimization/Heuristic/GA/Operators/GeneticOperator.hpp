#ifndef GENETICOPERATOR_HPP
#define GENETICOPERATOR_HPP

#include <vector>

namespace spgo
{
  template <typename Genotype>
  class GeneticOperator
  {
    public:
      GeneticOperator(){}
      virtual ~GeneticOperator(){}

      int getArity() { return arity; }

      virtual std::vector<Genotype> run(std::vector<Genotype>) = 0;

    protected:
      int arity;
  };
}

#endif