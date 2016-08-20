#ifndef REPLACEMENT_HPP
#define REPLACEMENT_HPP

#ifndef Gene_pool
#define Gene_pool std::vector<Genotype>
#endif

namespace spgo
{
  template<typename Fitness, typename Genotype, typename Phenotype>
  void SteadyState(
                   Gene_pool   population,
                   Gene_pool&  child     ,
                   Gene_pool   child_temp,
                   int*        selected  ,
                   double*     fitness   ,
                   int         arity     ,
                   Fitness     f         ,
                   Phenotype   (*grow)(Genotype)
                  )
  {
    for(int index = 0; index < arity; index++)
    {
      double evaluation;

      evaluation = f(grow(child_temp[index]));

      bool fitnes_flag = evaluation >= fitness[selected[index]];

      if(fitnes_flag)
      {
        child.push_back(child_temp[index]);
      }
      else
      {
        child.push_back(population[selected[index]]);
      }
    }
  }
}

#endif