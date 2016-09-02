#'@name simulate_dichotomous
#'@title Dichotomous data simulation
#'@param dim.data Data dimension
#'@param sample.size Size of the population that will be simulated
#'@param size.cluster Vector containing the number of items by dimension
#'@param seed_data Seed for simulation
#'@export 
simulate_dichotomous = function(dim.data = 1, sample.size = 1000, size.cluster = c(20), seed_data) {
  if ( dim.data > 1 )
    sim = simulate.poly.multi(sample.size = sample.size, size.cluster = size.cluster, dim.data = dim.data,
                              ncatgs = rep(2, sum(size.cluster)), seed_data = seed_data)  
  else 
    sim = simulate.poly.uni(n = sample.size, nitems = sum(size.cluster), 
                            ncatgs = rep(2, sum(size.cluster)), seed_data = seed_data)
  
  sim$data = sim$data - 1
  return (sim)
}
