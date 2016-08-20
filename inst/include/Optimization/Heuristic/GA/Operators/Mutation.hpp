#ifndef MUTATION_HPP
#define MUTATION_HPP

#include <vector>

namespace spgo
{
  template <typename Genotype>
  class Mutation : public GeneticOperator<Genotype>
  {
    public:
      Mutation() { this->arity = 1; }

      std::vector<Genotype> run(std::vector<Genotype> selected)
      {
        return mutate(selected[0]);
      }

    private:
      virtual std::vector<Genotype> mutate(Genotype) = 0;
  };
}

#endif