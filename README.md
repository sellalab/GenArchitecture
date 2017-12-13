# GenArchitecture
### This is the code for the polygenic trait simulations from Simons et al (2017).
### Running:
In this version, after 10N generations there is a shift in the optimum and we recored the response to the shift over 4000 generations. simulate.py runs the simulation, population.py includes classes for the population and mutations, and statWriter.py is a class for writing statistics to files.  
Thhis version recieves only three parameters (non optional) to be passed directly:  
shape of gamma distribution of scaled selection coefficients (default=1.0)  
scale of gamma distribution of scaled selection coefficients (default=10.0)  
shift in optimum 
  
Example: python simulate.py 10 2 0.2

## Branches
Branches provide variations on the main program. The main branch is the multidimensional version at equilibrium.
### OptimumShift - branch for testing the effects of a shift in the optimal phenotype (in 1D).
### MutBias - branch for testing the effects of mutational bias.

## Note:
Unlike the paper, the simulations use the convention that effect size is the difference between the two homozygotes.
The program was adapted to Python 3. A Python 2.7 version is available upon request.  
Many thanks to Guy Amster for helping write this code.
