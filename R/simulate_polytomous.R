#'@name simulate_polytomous
#'@title Polytomous data simulation
#'@param dim.data Data dimension
#'@param sample.size Size of the population that will be simulated
#'@param size.cluster Vector containing the number of items by dimension
#'@param ncatgs Vector containing the number of categories for each item
#'@param seed_data Seed for simulation
#'@export 
simulate_polytomous = function(dim.data = 1, sample.size = 1000, size.cluster = c(20), ncatgs = rep(4, 20), seed_data) {
  if ( dim.data > 1 )
    sim = simulate.poly.multi(sample.size = sample.size, size.cluster = size.cluster, dim.data = dim.data,
                              ncatgs = ncatgs, seed_data = seed_data)  
  else 
    sim = simulate.poly.uni(n = sample.size, nitems = sum(size.cluster), ncatgs = ncatgs, seed_data = seed_data)
  
  return (sim)
}
