from population import Mutation, MutationalProcess, Population
from collections import defaultdict, namedtuple
from statWriter import statWriter
from math import sqrt
import numpy as np
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
# frozenlist    :   list of mutations that were segregating before the shift, their effect size and their frequency before the shift and at present



# gather several statistics in one run
for shftFact in [0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0]:
    # set population with specific mutational and demographic parameters
    N       = 1000
    w       = 1.0
    u       = 1.0/100.0 
    shape   = float(sys.argv[1])
    scale   = float(sys.argv[2])
    mu = MutationalProcess(u, shape, scale)
    
    fitness = 'parents'
    pop = Population(N, w, mu, fitness)
    
    burnTime = 10*N 
    respTime = 4003 
    shift    = float(sys.argv[3])*shftFact
     
 
    ver = '2.2.6'
    # 2.2.6 list writing + phenotypic variation + lower mutation rate+N=1,000 again
    
    f = statWriter(os.getcwd()+'/results/',N=N,mu=mu,w=w,fitness=fitness,burnTime=burnTime,respTime=respTime,shift=shift,ver=ver)
    sc = 1.0 / float(2*N)
       
    sampleTimes = set([9000]+list(range(10000,10150,10))+list(range(10150,10500,50))+list(range(10500,12000,100))+list(range(12000,14000,250)))

    # we advance a total of 10*N+4000 generations 
    for time in range(burnTime+respTime):
        
        # advance one generation
        pop.nextGen()
        if time == burnTime:
            pop.shiftOptimum(shift)
            pop.freeze()
            i = 0
            for mut in pop.frozen():
                mut.index = i
                i+= 1
                #f.write('frozenList','init',mut.index,mut.scaledSize,mut.phenoSize, float(mut.frozenFreq) * sc)
        
        
        # once and then, we collect statistics
        if time in sampleTimes:
            print(time)

            meanF, stdF = pop.meanFitness()
            f.write('meanF',time,meanF)
            f.write('stdF',time,stdF)

            f.write('pheV',time,pop.pheVariance())
            
            seg = pop.segregating()
            f.write('numSeg',time,len(seg))
            
            f.write('zf',time,pop.zf)
            
            # collect statistics on segragating mutations
            meanZs = 0.0
            vsites=0.0
            f.write('segList','time',time)
            for mut in seg:
                p = float(mut.frequency) * sc
                q = 1.0 - p
                meanZs+= p*mut.phenoSize
                vsites+= 0.5*(mut.phenoSize**2)*p*q
                f.write('segList',mut.phenoSize, float(mut.frequency) * sc)
               
            
            f.write('meanZs',time,meanZs)
            f.write('stdZs',time,sqrt(vsites))

            if time>burnTime:
                frozen = pop.frozen()
                f.write('frozenList','time',time)
                for mut in frozen:
                    f.write('frozenList',mut.phenoSize, float(mut.frequency) * sc, float(mut.frozenFreq) * sc, mut.index)
                    

                    
                
                
            
    f.closeWriters()                


