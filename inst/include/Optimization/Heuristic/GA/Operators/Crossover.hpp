#ifndef CROSSOVER_HPP
#define CROSSOVER_HPP

#include <vector>

namespace spgo
{
  template <typename Genotype>
  class Crossover : public GeneticOperator<Genotype>
  {
    public:
      Crossover() { this->arity = 2; }

      std::vector<Genotype> run(std::vector<Genotype> selected)
      {
        return cross(selected[0], selected[1]);
      }

    private:
      virtual std::vector<Genotype> cross(Genotype, Genotype) = 0;
  };
}

#endif