# GenArchitecture
### This is the code for the polygenic trait simulations from Simons et al (2017).
### Running:
Program runs the full mltidimensional model and consists of three files. simulate.py runs the simulation, population.py includes classes for the population and mutations, and statWriter.py is a class for writing statistics to files.  
The program can recieve the following parameters:  
--n no. of traits (default=1)  
--N population size (default=1000)  
--w selection strength (default=1.0)  
--U mutation rate per haplotype genome (default=0.01)  
--shape shape of gamma distribution of scaled selection coefficients (default=1.0)  
--scale scale of gamma distribution of scaled selection coefficients (default=10.0)  
  
Example: python simulate.py --n 10 --U 0.1 --shape 10 --scale 2
__ __ __
## Branches
Branches provide variations on the main program. The main branch is the multidimensional version at equilibrium.
### OptimumShift - branch for testing the effects of a shift in the optimal phenotype (in 1D).
### MutBias - branch for testing the effects of mutational bias (in 1D).
__ __ __
## Submodules
Submodules provide additional files for the paper.
### Notebook - Provides a Mathematica notebook containing the main derivations and code to reproduce figures 1a, 2, 3 & 4.
### Inference - Provides Matlab files for fitting the high-pleiotropy, strong selection model to GWAS data. Includes data for height (from Wood et al., 2014) and BMI (from Locke et al., 2015). Used to produce figure 5.
### Demography - Provides C++ code to simulate the SFS under Schiffles & Durbin's MSMC-inferred demographic model. Used to produce figure 6.
__ __ __
## Note:
Unlike the paper, the simulations use the convention that effect size is the difference between the two homozygotes. 
The program was adapted to Python 3. A Python 2.7 version is available upon request.  
Many thanks to Guy Amster for helping write this code.
