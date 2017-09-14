# GenArchitecture
### This is the code for the polygenic trait simulations from Simons et al (2017).
### Running:
python simulate.py k Es n U  
Program expects 4 parameters - the shape and scale parameters of the distribution of scaled selection coefficients, the number of dimensions and the mutation rate per haploid genome.  
The program consists of three files. simulate.py runs the simulation, population.py includes classes for the population and mutations, and statWriter.py is a class for writing statistics to files.

## Branches
Branches provide variations on the main program. The main branch is the multidimensional version at equilibrium.
### Adaptation - branch for testing the effects of a shift in the optimal phenotype (in 1D).
### MutBias - branch for testing the effects of mutational bias (in 1D).

## Note:
Unlike the paper, the simulations use the convention the definition of effect size has the difference between the two homozygotes. 
