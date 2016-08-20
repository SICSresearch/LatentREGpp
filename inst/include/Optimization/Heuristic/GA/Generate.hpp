#ifndef GENERATE_HPP
#define GENERATE_HPP

#include <vector>

#ifndef Gene_pool
#define Gene_pool std::vector<Genotype>
#endif

namespace spgo
{
  template<typename Genotype>
  Gene_pool generateIndividuals(int pop_size)
  {
    Gene_pool result;
    for(int i = 0; i < pop_size; i++)
    {
      Genotype temp;
      result.push_back(temp);
    }

    return result;
  }
}

#endif