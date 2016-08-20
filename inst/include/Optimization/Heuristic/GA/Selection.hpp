#ifndef SELECTION_HPP
#define SELECTION_HPP

#include <random>

namespace spgo
{
  static void tournament(double* fitness, int size, int arity, int* result)
  {
    std::default_random_engine         generator(std::random_device{}());
    std::uniform_int_distribution<int> distribution(0, size - 1);

    int k = size/5;

    for(int j = 0; j < arity; j++)
    {
      int best = -1;

      for(int i = 0; i < k; i++)
      {
        int ind = distribution(generator);
        if(best == -1)
        {
          best = ind;
        }
        else
        {
          if(fitness[ind] > fitness[best])
          {
            best = ind;
          }
        }
      }
      result[j] = best;
    }
  }
}

#endif