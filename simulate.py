import argparse

from population import Mutation, MutationalProcess, Population
from collections import defaultdict, namedtuple
from statWriter import statWriter
from math import sqrt
import numpy as np
from scipy.stats import gamma, moment
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
# meanZ         :   mean phenotypic value
# meanZsites    :   mean expected genetic value - E(aq)
# zf            :   fixed contribution to phenotype
# stdZs         :   expected genetic variance - E(a^2 pq/2)
# skewZ         :   third moment of the phenotypic distribution

parser = argparse.ArgumentParser(description='Simulate population')
parser.add_argument('--N', type=int,help="population size",default=1000)
parser.add_argument('--w', type=float,help="selection strength",default=1.0)
parser.add_argument('--U', type=float,help="mutation rate",default=0.01)
parser.add_argument('--shape', type=float,help="shape",default=1.0)
parser.add_argument('--biases', type=float,help="biases",default=[0], nargs='+')
parser.add_argument('--scale', type=float,help="scale",default=10)
parser.add_argument('--pweak', type=float,help="fraction weakly selected",default=0.5)

args=parser.parse_args()

print(args)


# gather several statistics in one run
for mutFact in args.biases:
    # set population with specific mutational and demographic parameters
    N       = args.N
    w       = args.w
    u       = args.U
    n       = 1 #Only 1D. MultiD not implemented.
    shape   = args.shape
    scale   = args.scale
    mu = MutationalProcess(u, shape, scale)
    bias=mutFact
    
    fitness = 'parents'
    pop = Population(N, w, mu, fitness, args.pweak)
    
    burnTime = 10*N 
    respTime = 10003 
    shift    = 0
   
     
 
    ver = '2.2.14'
    # 2.2.14  Mutation pressu
    #/Users/ybs2103/PhD/Python/Adaptation/Results/
    #/ifs/data/c2b2/gs_lab/shared/ybs2103/Adaptation/results/

    f = statWriter(os.getcwd()+'/results',N=N,mu=mu,w=w,fitness=fitness,burnTime=burnTime,respTime=respTime,shift=shift,ver=ver,bias=bias,pweak=args.pweak)
    sc = 1.0 / float(2*N)
       
    sampleTimes = set(list(range(0,10000,1000))+list(range(10000,20001,1000)))

    # we advance a total of 10*N+10000 generations 
    for time in range(burnTime+respTime):
        
        if (time % 100)==0:
            print(time)

        # advance one generation
        pop.nextGen()
        if time == burnTime:
            pop.shiftOptimum(shift)
            pop.setBias(bias)
        
        
        # once and then, we collect statistics
        if time in sampleTimes:
            print(time)
            f.write('test',1.0,2.0)

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
            phenos=pop.phenotypes()
            if time==20000:
                #f.write('segList','time',time)
                #f.write('pheList','time',time)
                
                for phe in  phenos:
                    f.write('pheList',phe)

                #for pheno in pop.denovo:
                #    f.write('denovo',pheno) 
            for mut in seg:
                p = float(mut.frequency) * sc
                q = 1.0 - p
                meanZs+= p*mut.phenoSize
                vsites+= 0.5*(mut.phenoSize**2)*p*q
                if time==20000:
                    f.write('segList',mut.phenoSize, float(mut.frequency) * sc,mut.ES)
                    

            f.write('meanZ',time,np.mean(phenos))
            f.write('meanZsites',time,pop.zf+meanZs)
            f.write('Zf',time,pop.zf)
            f.write('stdZs',time,sqrt(vsites))
            f.write('skewZ',time,moment(phenos,3))

                    
                
                
            
    f.closeWriters()                


