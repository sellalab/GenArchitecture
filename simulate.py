import time as zman
t1=zman.time()
import argparse
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


parser = argparse.ArgumentParser(description='Simulate population')
parser.add_argument('--n', type=int,help="no. of traits",default=1)
parser.add_argument('--N', type=int,help="population size",default=1000)
parser.add_argument('--w', type=float,help="selection strength",default=1.0)
parser.add_argument('--U', type=float,help="mutation rate",default=0.01)
parser.add_argument('--shape', type=float,help="shape",default=1.0)
parser.add_argument('--scale', type=float,help="scale",default=10)


args=parser.parse_args()

print(args)

# set population with specific mutational and demographic parameters
N       = args.N
w       = args.w
u       = args.U
n       = args.n  
shape   = args.shape
scale   = args.scale
mu = MutationalProcess(u, shape, scale)


fitness = 'parents'
pop = Population(N, w, mu,n, fitness)

burnTime = 10*N 
respTime = 3 
shift    = 0
   
 
 
ver = 'mD 1.0'
# mD 1.0 multidimensional
#/Users/ybs2103/PhD/Python/Adaptation/Results/
#/ifs/data/c2b2/gs_lab/shared/ybs2103/Adaptation/results/

f = statWriter(os.getcwd()+'/results'+str(N),N=N,mu=mu,n=n,w=w,fitness=fitness,burnTime=burnTime,respTime=respTime,shift=shift,ver=ver)
sc = 1.0 / float(2*N)

print(os.getcwd()+'/results')

sampleTimes = set(list(range(0,11*N,N//2)))

# we advance a total of 20*N generations 
for time in range(burnTime+respTime):
    
    if (time % (N/10))==0:
        print("Generation no. " + str(time) + " is " + str(time/(N/10.0)) + "%")

    # advance one generation
    pop.nextGen()    
    
    # once and then, we collect statistics
    if time in sampleTimes:
        print(time)

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
        f.write('meanZ',time,pop.zf+meanZs)
        
        f.write('realTime',time, (zman.time()-t1)/3600)
        
        f.flush()

                
            
            
        
f.closeWriters()                


