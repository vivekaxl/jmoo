
"""
##########################################################
### @Author Joe Krall      ###############################
### @copyright see below   ###############################

    This file is part of JMOO,
    Copyright Joe Krall, 2014.

    JMOO is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    JMOO is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with JMOO.  If not, see <http://www.gnu.org/licenses/>.
    
###                        ###############################
##########################################################
"""

"Brief notes"
"Standardized MOEA code for running any MOEA"

from jmoo_algorithms import *
from jmoo_stats_box import *
from jmoo_properties import *
from Moo import *
from pylab import *

def jmoo_evo(problem, algorithm, toStop = bstop):
    """
    ----------------------------------------------------------------------------
     Inputs:
      -@problem:    a MOP to optimize
      -@algorithm:  the MOEA used to optimize the problem
      -@toStop:     stopping criteria method
    ----------------------------------------------------------------------------
     Summary:
      - Evolve a population for a problem using some algorithm.
      - Return the best generation of that evolution
    ----------------------------------------------------------------------------
     Outputs:
      - A record (statBox) of the best generation of evolution
    ----------------------------------------------------------------------------
    """
    
    # # # # # # # # # # #
    # 1) Initialization #
    # # # # # # # # # # #
    stoppingCriteria = False                             # Just a flag for stopping criteria
    statBox          = jmoo_stats_box(problem,algorithm) # Record keeping device
    gen              = 0                                 # Just a number to track generations
    
    # # # # # # # # # # # # # # # #
    # 2) Load Initial Population  #
    # # # # # # # # # # # # # # # #
    population = problem.loadInitialPopulation(MU)
    #Vivek:
    #print "Length of population: ",len(population)
    
    # # # # # # # # # # # # # # #
    # 3) Collect Initial Stats  #
    # # # # # # # # # # # # # # #
    statBox.update(population, 0, 0, initial=True)
    
    # # # # # # # # # # # # # # #
    # 4) Generational Evolution #
    # # # # # # # # # # # # # # #
    
    while gen < PSI and stoppingCriteria == False:
        gen+= 1
        #print "Generation: ",gen
        # # # # # # # # #
        # 4a) Selection #
        # # # # # # # # #
        if algorithm.name == "GALE2_1" and gen == 1:
            # use the initial population to build classes
            for pop in population:
                pop.fitness = problem.evaluate(pop.decisionValues)
            high = -1e6
            constraints = []
            for x in xrange(len(problem.decisions)):
                for e,d in enumerate(sdiv2(population, cohen=0.3, num1=lambda y: y.decisionValues[x], num2=lambda y: min(y.fitness))):
                    temp = sorted([y.decisionValues[x] for y in d[-1]])
                    mean =  sum([min(y.fitness) for y in d[-1]])/len(d[-1])
                    if mean > high:
                        const1 = temp[0]
                        const2 = temp[-1]
                        high = mean
                problem.decisions[x].low = const1
                problem.decisions[x].up = const2
            # read the data from population
            # run sdiv get constraints and input it in the problem. There is a problem with upper limit

            
        problem.referencePoint = statBox.referencePoint
        selectees,evals = algorithm.selector(problem, population)
        numNewEvals = evals
        
        
        #raw_input("Press any Key")
        # # # # # # # # # #
        # 4b) Adjustment  #
        # # # # # # # # # #
        selectees,evals = algorithm.adjustor(problem, selectees)
        numNewEvals += evals
        
        # # # # # # # # # # #
        # 4c) Recombination #
        # # # # # # # # # # #
        population,evals = algorithm.recombiner(problem, population, selectees, MU)
        numNewEvals += evals        
        
        
        # # # # # # # # # # #
        # 4d) Collect Stats #
        # # # # # # # # # # #
        statBox.update(population, gen, numNewEvals)
        #for row in population: print row.valid
        
        
            
        # # # # # # # # # # # # # # # # # #
        # 4e) Evaluate Stopping Criteria  #
        # # # # # # # # # # # # # # # # # #
        stoppingCriteria = toStop(statBox)
        #print ">>>>>>>> Stopping Criteria: ",stoppingCriteria
        #stoppingCriteria = False
        
        
    


    #return the representative generation
    return statBox


def sdiv2(lst, tiny=3, cohen=0.3, num1=lambda x: x[0], num2=lambda x: x[1]):
        "Divide lst of (num1,num2) using variance of num2."
        #----------------------------------------------
        class Counts():  # Add/delete counts of numbers.
            def __init__(i, inits=[]):
                i.zero()
                for number in inits: i + number

            def zero(i):
                i.n = i.mu = i.m2 = 0.0

            def sd(i):
                if i.n < 2:
                    return i.mu
                else:
                    return (max(0, i.m2) * 1.0 / (i.n - 1)) ** 0.5

            def __add__(i, x):
                i.n += 1
                delta = x - i.mu
                i.mu += delta / (1.0 * i.n)
                i.m2 += delta * (x - i.mu)

            def __sub__(i, x):
                if i.n < 2: return i.zero()
                i.n -= 1
                delta = x - i.mu
                i.mu -= delta / (1.0 * i.n)
                i.m2 -= delta * (x - i.mu)
                #----------------------------------------------

        def divide(this, small):  #Find best divide of 'this'
            lhs, rhs = Counts(), Counts(num2(x) for x in this)
            n0, least, cut = 1.0 * rhs.n, rhs.sd(), None
            for j, x in enumerate(this):
                if lhs.n > tiny and rhs.n > tiny:
                    maybe = lhs.n / n0 * lhs.sd() + rhs.n / n0 * rhs.sd()
                    if maybe < least:
                        if abs(lhs.mu - rhs.mu) >= small:  # where's the paper for this method?
                            cut, least = j, maybe
                rhs - num2(x)
                lhs + num2(x)
            return cut, least

        #----------------------------------------------
        def recurse(this, small, cuts):
            #print this,small
            cut, sd = divide(this, small)
            if cut:
                recurse(this[:cut], small, cuts)
                recurse(this[cut:], small, cuts)
            else:
                cuts += [(sd, this)]
            return cuts

        #---| main |-----------------------------------
        # for x in lst:
        #   print num2(x)
        small = Counts(num2(x) for x in lst).sd() * cohen  # why we use a cohen??? how to choose cohen
        if lst:
            return recurse(sorted(lst, key=num1), small, [])

