
Code used to run algorithm generate final figures. 
Algorithm first generates random phylogenic tree. Raw or “observed” data is then generated with error rates using Poisson and multinomial distributions. Using Bayesian statistic, a posterior probability is calculated of the genotype at each position at each location. Mutation statues are assigned to each location and clustering is performed for each cell to build a phylogenic tree of the cells’ evolution. With each mutation we calculated the posterior probability of where it most likely occurred on the evolutionary tree.
