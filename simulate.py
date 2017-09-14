from population import Mutation, MutationalProcess, Population
from collections import defaultdict, namedtuple
from statWriter import statWriter
from math import sqrt
import numpy as np
np.set_printoptions(threshold='nan',linewidth=2500000)
from scipy.stats import gamma
import sys
import os

# statistical summary of the population state in a given point of time
# time          :   generation number.
# meanF         :   mean fitness of the population.
# stdF          :   stadrard deviation of population fitness.
# pheV          :   the pheotypic variance.
# numSeg        :   number of segregating mutations.
# segList       :   List of effect sizes and frequencies of segregating mutations.
# pheList       :   List of phenotypic values of individuals in the populations.
# meanZs        :   mean expected genetic value - E(aq)
# stdZs         :   expected genetic variance - E(a^2 pq/2)




# gather several statistics in one run
for mutFact in [0.0]:
    # set population with specific mutational and demographic parameters
    N       = 1000
    w       = 1.0
    shape   = float(sys.argv[1])
    scale   = float(sys.argv[2])
    n=int(sys.argv[3])
    u       = float(sys.argv[4]) 

    mu = MutationalProcess(u, shape, scale)
    bias=0.0

    
    fitness = 'parents'
    pop = Population(N, w, mu,n, fitness)
    
    burnTime = 10*N 
    respTime = 3 
    shift    = 0
   
     
 
    ver = 'mD 1.0'
    # mD 1.0 multidimensional
    #/Users/ybs2103/PhD/Python/Adaptation/Results/
    #/ifs/data/c2b2/gs_lab/shared/ybs2103/Adaptation/results/

    f = statWriter(os.getcwd()+'/results',N=N,mu=mu,n=n,w=w,fitness=fitness,burnTime=burnTime,respTime=respTime,shift=shift,ver=ver,bias=bias)
    sc = 1.0 / float(2*N)
    
    print os.getcwd()+'/results'

    sampleTimes = set(range(0,11*N,N/2))

    # we advance a total of 10*N+3 generations 
    for time in xrange(burnTime+respTime):
        
        if (time % (N/100))==0:
            print("Generation no. " + str(time) + " is " + str(time/(N/10.0)) + "%")

        # advance one generation
        pop.nextGen()
        if time == burnTime:
            pop.shiftOptimum(shift)
            pop.setBias(bias)
        
        
        # once and then, we collect statistics
        if time in sampleTimes:
            print time

            meanF, stdF = pop.meanFitness()
            f.write('meanF',time,meanF)
            f.write('stdF',time,stdF)

            f.write('pheV',time,pop.pheVariance())
            
            seg = pop.segregating()
            f.write('numSeg',time,len(seg))
            
            #f.write('zf',time,pop.zf)
            
            # collect statistics on segragating mutations
            meanZs = 0.0
            vsites=0.0
            if time==10*N:
                f.write('segList','time',time)
                #f.write('pheList','time',time)
                phenos=pop.phenotypes()
                for phe in  phenos:
                    f.write('pheList',phe)

                #for pheno in pop.denovo:
                #    f.write('denovo',pheno) 
            for mut in seg:
                p = float(mut.frequency) * sc
                q = 1.0 - p
                meanZs+= p*mut.phenoSize
                vsites+= 0.5*(mut.phenoSize**2)*p*q
                if time==10*N:
                    f.write('segList',mut.phenoSize, float(mut.frequency) * sc)
                    

            f.write('meanZ',time,pop.zf+meanZs)
            f.write('stdZs',time,np.sqrt(vsites))


                    
                
                
            
    f.closeWriters()                


