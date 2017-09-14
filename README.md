# GenArchitecture
### This is the code for the polygenic trait simulations from Simons et al (2017).
### Running:
This version has multiple parameters and so uses Python’s argument parser.  
Example: simulate.py —-U 3 —-shape 10 —biases 0 2 —scale 3 —pweak 0.5  
Though it’s included as a parameter for future use, the number of traits has to be left as 1 since the multidimensional version is yet to be implemented.  
The program consists of three files. simulate.py runs the simulation, population.py includes classes for the population and mutations, and statWriter.py is a class for writing statistics to files.

## Branches
Branches provide variations on the main program. The main branch is the multidimensional version at equilibrium.
### Adaptation - branch for testing the effects of a shift in the optimal phenotype (in 1D).
### MutBias - branch for testing the effects of mutational bias (in 1D).

## Notes:
Unlike the paper, the simulations use the convention the definition of effect size has the difference between the two homozygotes.  
Code has been adapted to python 3. A python 2.7 version is available upon request.