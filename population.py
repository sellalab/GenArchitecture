from collections import defaultdict, namedtuple
import math
import random
import numpy as np

# Container class for mutations.
class Mutation(object):
    def __init__(self, scaledSize, phenoSize):
        self.scaledSize = scaledSize
        self.phenoSize  = phenoSize
        self.frequency  = 1
        self.frozenFreq = -1 # start undefined
        self.index      = -1 # start undefined
            

# parameters for the mutational process
# mu: haploid mutation rate per-generation
# shape, scale: gamma-distribution parameters, describing the effect size distribution
MutationalProcess = namedtuple('MutationalProcess', 'mu shape scale')

class Population(object):
    
    #   N : population size
    def __init__(self, N, w, mu, n, selectionMode):
        
        # population size
        self.N = N

        # dimensionality
        self.n = n        
        
        # set of segregating mutations
        self._segregating = set()

        # set of de-novo mutations
        self.denovo = set()
        
        # list of mutations carried by each individual of the population:
        # self._individuals[i] is a dictionary representing the genotype of individual i in the population:
        #   - The keys are mutations present in the individual
        #   - Values are the ploidity of the mutation (either 1 or 2)
        self._individuals = [defaultdict(int) for _ in range(self.N)]
        
        # same data structure as self._individuals is used as a temporary data-structure when constructing next generation
        self._offspring = [defaultdict(int) for _ in range(self.N)]
                
        # total fitness effect of fixed mutations in the population
        self.zf = np.zeros(n)
        
        # initialize fitness optimum as zero
        self._fitnessOptimum = np.zeros(n)
        
        # selection model paramter
        self._fitnesscoef = 0.5/float(w*w)
        
        # mutation scaling parameter
        self._muScalingCoef = 2.0*float(w*w)/float(N)
        
        # mutational process parameters
        self._mu = mu
        
        # segregation update counter
        self._lastUpdate = 0

        # mutation bias
        self._bias=0
        
        # read selection mode
        if selectionMode == 'parents':
            self._sampleOffspring = self._sampleOffspringByParentsFitness
        elif selectionMode == 'offspring':
            self._sampleOffspring = self._sampleOffspringByChildrenFitness
            
    def shiftOptimum(self, delta):
        self._fitnessOptimum+= delta

    def setBias(self, bias):
        self._bias= bias
    
    def nextGen(self):
        
        # increment generation counter
        self._lastUpdate+= 1
        
        # create all offspring:
        self._sampleOffspring()        
        # advance generation:
        
        # swap individuals and offspring
        self._individuals, self._offspring = self._offspring, self._individuals
            
        # once every few iterations, update the segregating mutations list:
        if self._lastUpdate >= 64:
            self._updateSegregatingList()
    
    def segregating(self):
        # update the segregating mutations list
        self._updateSegregatingList()
        
        return self._segregating
    
    def freeze(self):
        self._frozen = set()
        for mu in self.segregating():
            self._frozen.add(mu)
            mu.frozenFreq = mu.frequency
    
    def frozen(self):
        self._updateSegregatingList()
        
        return self._frozen
    
    # This method updates the segregating mutations list:
    # - Compute the frequency of each mutation
    # - Removes extinct or fixed mutations
    
    def _updateSegregatingList(self):
        
        if self._lastUpdate == 0:
            return
        
        # reset counts
        for mu in self._segregating:
            mu.frequency = 0
        
        # update counts
        for indv in self._individuals:
            for mu,ploidity in iter(indv.items()):
                mu.frequency+= ploidity
                
        # iterate over segregating mutation list and search for extinct or fixed mutations
        extinct = [mu for mu in self._segregating if mu.frequency==0]
        fixed   = [mu for mu in self._segregating if mu.frequency==(2*self.N)]
        
        # for fixed mutations
        for mu in fixed:
            
            # remove mutation from all offspring:
            for indv in self._individuals:
                del indv[mu]

            # add mutation effect size to total fixed effect
            self.zf += mu.phenoSize

        # remove extinct and fixed mutations from the segregating list:
        for mu in (extinct + fixed):
            self._segregating.remove(mu)
        
        # reset last-update time
        self._lastUpdate = 0
    
    
    # compute the fitness of an individual
    # indv : input individual (segregating genotype dictionary)
    def _fitness(self, indv):
        
        z = self._phenotype(indv)  - self._fitnessOptimum
        
        return math.exp(-np.dot(z,z)*self._fitnesscoef)


    # get phenotype
    def _phenotype(self, indv):
        
        z = 0.0
        
        # sum the effect of segregating mutations 
        for mu,ploidity in iter(indv.items()):
            z+= (ploidity*mu.phenoSize)
            
        # add the effect of fixed mutations, and substract the optimal value
        z = self.zf + 0.5*z
        
        return z
    
    def phenotypes(self):
        return [self._phenotype(self._individuals[i]) for i in range(self.N)]

    # sample self._offspring, based on parents fitness
    def _sampleOffspringByParentsFitness(self):
        
        # we will sample parents with probabilities proportional to their fitness.
        # first, we compute the fitness
        fitVals = [self._fitness(self._individuals[i]) for i in range(self.N)]
        # we normalize to a probability vector and compute the CDF.
        tmp = 1.0 / np.sum(fitVals)
        cdf = np.add.accumulate([x*tmp for x in fitVals])
        
        # generate offspring one by one
        for child in self._offspring:
            
            # sample two non-identical random parents from the cdf
            p1, p2 = np.searchsorted(cdf,random.random()), np.searchsorted(cdf,random.random())
            while p1 == p2:
                p2 = np.searchsorted(cdf,random.random())
            
            # generate random child of these two parents
            self._sampleChild(child, self._individuals[p1], self._individuals[p2])
    
    # return the mean fitness oand std f the population
    def meanFitness(self):
        fit = np.array([self._fitness(self._individuals[i]) for i in range(self.N)])
        return np.mean(fit), np.std(fit) 

    # phenotypic variance
    def pheVariance(self):
        phe = np.array([self._phenotype(self._individuals[i]) for i in range(self.N)])
        return np.var(phe,axis=0)


  
     # sample self._offspring, based on offspring fitness
    def _sampleOffspringByChildrenFitness(self):
        
        # generate offspring one by one
        index = 0
        while index < self.N:
                        
            # sample distinct random parents uniformly
            p1, p2 = random.randint(0,self.N - 1), random.randint(0,self.N - 1)
            
            while p1 == p2:
                p2 = random.randint(0,self.N - 1)
            
            # generate random child of these two parents
            self._sampleChild(self._offspring[index], self._individuals[p1], self._individuals[p2])
            
            # reject or accept the resulting child with probability proportional to its fitness
            if random.random() < self._fitness(self._offspring[index]) :
                index+= 1

    
    # samples one random child of parents p1 and p2
    def _sampleChild(self, child, p1, p2):
        
        # clear offspring dictionary
        child.clear()
        
        # for each parent
        for p in [p1, p2]:
            
            # for each mutation carried by the parent
            for mu,ploidity in iter(p.items()):
                
                # if the parent has two copies of the mutation he is bound to pass it on
                if ploidity == 2:
                    child[mu] += 1
                
                # if the parent is heterozygous, he pass it on with probabilty 0.5
                elif random.getrandbits(1):
                    child[mu] += 1
        
        #self.denovo.clear()

        # add random de-novo mutations:
        # number of de-novo mutations is a poisson variable
        for _ in range(np.random.poisson(2.0*self._mu.mu)):
            
            # the scaled effect size has gamma distribution
            scaledSize = np.random.gamma(self._mu.shape, self._mu.scale)
            v=np.random.multivariate_normal(np.zeros(self.n),np.identity(self.n))
            phenoSize  = v*math.sqrt(self._muScalingCoef*scaledSize/np.dot(v,v))
            

            #self.denovo.add(phenoSize)

            # we add the mutation to the segregating list and to the new offspring (in heterozygous state)
            mu = Mutation(scaledSize,phenoSize)
            self._segregating.add(mu)
            child[mu] = 1
    

