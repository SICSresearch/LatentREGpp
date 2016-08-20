#ifndef GA_HPP
#define GA_HPP

#include <Optimization/Optim.hpp>
#include <Optimization/Heuristic/GA/Operators/GeneticOperator.hpp>
#include <Optimization/Heuristic/GA/Generate.hpp>
#include <cmath>
#include <random>
#include <vector>
/////////////////////////
#include <iostream>
/////////////////////////

#define Gene_pool   std::vector<Genotype>
#define Generate    Gene_pool(*generate_function)(int)
#define Grow        Phenotype(*grow_function)(Genotype)
#define Selection   void(*selection_function)(double*, int, int, int*)
#define Replacement void(*replacement_function)(Gene_pool, Gene_pool&, Gene_pool, int*, double*, int, Fitness, Grow)

namespace spgo
{
  template <typename Fitness, typename Genotype, typename Phenotype>
  class GeneticAlgorithm : public Optim<Fitness, Phenotype, void*>
  {
    public:

      GeneticAlgorithm(Generate, Grow, Selection, Replacement, int pop_size, std::vector<GeneticOperator<Genotype>*> operators)
      {
        init(generate_function, grow_function, selection_function, pop_size, 10000, operators);
      }

      GeneticAlgorithm(Grow,
                       Selection,
                       Replacement,
                       int pop_size,
                       int generations_limit,
                       std::vector<GeneticOperator<Genotype>*> operators)
      {
        init(generateIndividuals,
             grow_function,
             selection_function,
             replacement_function,
             pop_size,
             generations_limit,
             operators);
      }

      void init(Generate,
                Grow,
                Selection,
                Replacement,
                int pop_size,
                int generations_limit,
                std::vector<GeneticOperator<Genotype>*> operators)
      {
        this->generate_function    = generate_function;
        this->grow_function        = grow_function;
        this->selection_function   = selection_function;
        this->replacement_function = replacement_function;
        this->pop_size             = pop_size;
        this->generations_limit    =generations_limit;
        fitness_eval               = new double[pop_size];

        Operators.push_back(operators[0]);
        Operators.push_back(operators[1]);
      }

      ~GeneticAlgorithm()
      {
        delete[] fitness_eval;
      }

      int run(std::vector<Fitness> f, Phenotype &p, void*)
      {
        std::vector<int>                   index;
        std::random_device                 rd;
        std::mt19937                       e2(rd());
        std::uniform_int_distribution<int> dist(0,100);

        // Initialization
        population = (*generate_function)(pop_size);

        // Evaluation
        double max = -9999;
        double min =  9999;
        double sum =  0;

        for(int i = 0; i < pop_size; i++)
        {
          fitness_eval[i] = (f[0])(grow_function(population[i]));
          if(fitness_eval[i] > max) { max = fitness_eval[i]; }
          if(fitness_eval[i] < min) { min = fitness_eval[i]; }
          sum += fitness_eval[i];
        }

        std::cout << "-1 | " << max << " | " << min << " | " << sum/pop_size << std::endl;

        for(int i = 0; i < generations_limit; i++)
        {
          bool      finish_flag;
          Gene_pool child;

          finish_flag = false;
          max         = -9999;
          min         =  9999;
          sum         =  0;

          while(!finish_flag)
          {
            int       gen_op;
            int*      index;
            Gene_pool selected;

            gen_op = dist(e2) >= 30.0 ? 0 : 1;

            index  = new int[Operators[gen_op]->getArity()];

            selection_function(fitness_eval,
                               pop_size,
                               Operators[gen_op]->getArity(),
                               index);

            finish_flag = child.size() >= pop_size;

            for(int j = 0; j < Operators[gen_op]->getArity(); j++)
            {
              selected.push_back(population[index[j]]);
            }

            Gene_pool child_temp = Operators[gen_op]->run(selected);

            replacement_function(population,
                                 child,
                                 child_temp,
                                 index,
                                 fitness_eval,
                                 Operators[gen_op]->getArity(),
                                 f[0],
                                 grow_function);

            delete[] index;
          }

          for(int idx = 0; idx < pop_size; idx++)
          {
            population[idx] = child[idx];
            fitness_eval[idx] = (f[0])(grow_function(population[idx]));
            if(fitness_eval[idx] > max) { max = fitness_eval[idx]; }
            if(fitness_eval[idx] < min) { min = fitness_eval[idx]; }
            sum += fitness_eval[idx];
          }

          std::cout << i << " | " << max << " | " << min << " | " << sum/pop_size << std::endl;
        }

        max         = -9999;
        min         =  9999;
        sum         =  0;

        for(int i = 0; i < pop_size; i++)
        {
          fitness_eval[i] = (f[0])(grow_function(population[i]));
          if(fitness_eval[i] > max) { max = fitness_eval[i]; }
          if(fitness_eval[i] < min) { min = fitness_eval[i]; }
          sum += fitness_eval[i];
        }

        min = 99999;

        for(int i = 0; i < pop_size; i++)
        {
          Phenotype temp_p = grow_function(population[i]);
          double temp = (f[0])(temp_p);

          if(temp < min)
          {
            min = temp;
            p   = temp_p;
          }
        }

        std::cout << "final | " << max << " | " << min << " | " << sum/pop_size << std::endl;

        return 0;
      }

      private:
      	Generate;
        Grow;
        Gene_pool   population;
        Selection;
        Replacement;

        double*     fitness_eval;
        int         pop_size,
                    generations_limit;

        std::vector<GeneticOperator<Genotype>*> Operators;
    };
}

#endif