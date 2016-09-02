#'@name simulate_polytomous
#'@title data simulation
#'@export
simulate_polytomous = function(dim.data = 1, sample.size = 1000, size.cluster = c(20), ncatgs = rep(4, 20), seed_data) {
  if ( dim.data > 1 )
    sim = simulate.poly.multi(sample.size = sample.size, size.cluster = size.cluster, dim.data = dim.data,
                              ncatgs = ncatgs, seed_data = seed_data)  
  else 
    sim = simulate.poly.uni(n = sample.size, nitems = sum(size.cluster), ncatgs = ncatgs, seed_data = seed_data)
  
  return (sim)
}
