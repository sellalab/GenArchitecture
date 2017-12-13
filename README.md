# GenArchitecture
### This is the code for the polygenic trait simulations from Simons et al (2017).
### Running:
Version for running analysis with mutational bias. Here, new mutations are either strongly selected with E(S)=50 or weakly selected with E(S)=0.1. Therefore, the scale parameter is ignored and we introduce a parameter for the proportion of weakly selected new mutations.  simulate.py runs the simulation, population.py includes classes for the population and mutations, and statWriter.py is a class for writing statistics to files.  
The program can recieve the following parameters:    
--N population size (default=1000)  
--w selection strength (default=1.0)  
--U mutation rate per haplotype genome (default=0.01)  
--shape shape of gamma distribution of scaled selection coefficients (default=1.0)  
[ignored] --scale scale of gamma distribution of scaled selection coefficients (default=10.0)  
--biases list of mutation biases to run (default=[0])  
--pweak proportion of new mutations that are weakly selected   
  
Example: python simulate.py --n 10 --U 0.1 --biases 0 0.2 0.3 --pweak 0.3

## Branches
Branches provide variations on the main program. The main branch is the multidimensional version at equilibrium.
### OptimumShift - branch for testing the effects of a shift in the optimal phenotype (in 1D).
### MutBias - branch for testing the effects of mutational bias (in 1D).

## Note:
Unlike the paper, the simulations use the convention that effect size is the difference between the two homozygotes.
The program was adapted to Python 3. A Python 2.7 version is available upon request.  
Many thanks to Guy Amster for helping write this code.
