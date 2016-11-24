#'@name simulate_dichotomous
#'@title Dichotomous data simulation 
#'@description This function generates dichotomous test 
#'@param dim.data Data dimension
#'@param sample.size Size of the population to be simulated
#'@param size.cluster Vector containing the number of items per dimension
#'@param seed_data Seed for simulation of the data
#'@param seed_item Seed for simulation of the item parameters
#'@param model the model used to generate the test. It can have values of "1PL", "2PL" or "3PL".
#'@return list with the simulated data, the information of the item parameters and
#'the individual parameters
#'@examples
#'
#'#Example 1: see by default options
#'simm=simulate_dichotomous()
#'@export 

simulate_dichotomous = function(dim.data = 1, sample.size = 1000, size.cluster = c(20), seed_data=5000L,seed_item=1000,model="2PL") {
  if ( dim.data > 1 )
    sim = simulate.poly.multi(sample.size = sample.size, size.cluster = size.cluster, dim.data = dim.data,
                              ncatgs = rep(2, sum(size.cluster)), seed_data = seed_data)  
  else 
    sim = simulate.poly.uni(n = sample.size, nitems = sum(size.cluster), 
                            ncatgs = rep(2, sum(size.cluster)), seed_data = seed_data,
                            seed_item=seed_item,model=model)
  
  sim$data = sim$data - 1
  return (sim)
}
