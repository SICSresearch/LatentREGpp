#'@name simulate_polytomous
#'@title Polytomous data simulation
#'@description This function generates dichotomous test of an specified dimension and size. The items per dimension are specified in the size.cluster vector, and the individuals are specified in sample.size.
#'The amout of categories of each item is determined with ncatgs.
#'@param dim.data Data dimension
#'@param sample.size Size of the population that will be simulated
#'@param size.cluster Vector containing the number of items by dimension
#'@param ncatgs Vector containing the number of categories for each item
#'@param seed_data Seed for simulation
#'@export 
simulate_polytomous = function(dim.data = 1, sample.size = 1000, size.cluster = c(20), ncatgs = rep(4, 20), seed_data=5000L,model="2PL") {
  if ( dim.data > 1 )
    sim = simulate.poly.multi(sample.size = sample.size, size.cluster = size.cluster, dim.data = dim.data,
                              ncatgs = ncatgs, seed_data = seed_data)  
  else 
    sim = simulate.poly.uni(n = sample.size, nitems = sum(size.cluster), ncatgs = ncatgs, seed_data = seed_data,model=model)
  
  return (sim)
}
