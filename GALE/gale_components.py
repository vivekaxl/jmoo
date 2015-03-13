"""
    This file is part of GALE,
    Copyright Joe Krall, 2014.

    GALE is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    GALE is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with GALE.  If not, see <http://www.gnu.org/licenses/>.
"""

import os, sys, inspect, random

cmd_subfolder = os.path.realpath(os.path.abspath(os.path.join(os.path.split(inspect.getfile( inspect.currentframe()))[0],"fastmap")))
if cmd_subfolder not in sys.path:
    sys.path.insert(0, cmd_subfolder)
    
from Slurp import *
from Moo import *
from Moo2 import *
from jmoo_individual import *

currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir)

import jmoo_properties

def galeWHERE(problem, population):
    "The Core method behind GALE"
    
    # Compile population into table form used by WHERE
    # t is a table which has rows and each row has cells
    t = slurp([[x for x in row.decisionValues] + ["?" for y in problem.objectives] for row in population], problem.buildHeader().split(","))


    # Initialize some parameters for WHERE
    The.allowDomination = True
    The.alpha = 1
    for i,row in enumerate(t.rows): # Attention JOE: at this point nothing is evaluated, but there might be some rows which are evaluated
        row.evaluated = False
    
    # Run WHERE
    m = Moo(problem, t, len(t.rows), N=1)
    m.divide(minnie=rstop(t))
          
    # Organizing
    NDLeafs = m.nonPrunedLeaves()                       # The surviving non-dominated leafs
    allLeafs = m.nonPrunedLeaves() + m.prunedLeaves()   # All of the leafs
    
    # After mutation: Check how many rows were actually evaluated
    numEval = 0
    for leaf in allLeafs:
        for row in leaf.table.rows:
            if row.evaluated:
                numEval += 1
                
    return NDLeafs, numEval

def galeWHERE2(problem, population):
    "The Core method behind GALE"

    # Compile population into table form used by WHERE
    # t is a table which has rows and each row has cells
    t = slurp([[x for x in row.decisionValues] + ["?" for y in problem.objectives] for row in population], problem.buildHeader().split(","))


    # Initialize some parameters for WHERE
    The.allowDomination = True
    The.alpha = 1
    for i,row in enumerate(t.rows): # Attention JOE: at this point nothing is evaluated, but there might be some rows which are evaluated
        row.evaluated = False

    # Run WHERE
    m = Moo2(problem, t, len(t.rows), N=1)
    m.divide(minnie=rstop(t))

    # Organizing
    NDLeafs = m.nonPrunedLeaves()                       # The surviving non-dominated leafs
    allLeafs = m.nonPrunedLeaves() + m.prunedLeaves()   # All of the leafs

    # After mutation: Check how many rows were actually evaluated
    numEval = 0
    for leaf in allLeafs:
        for row in leaf.table.rows:
            if row.evaluated:
                numEval += 1

    return NDLeafs, numEval
      
     
def galeMutate(problem, NDLeafs):
    
    #################
    # Mutation Phase
    #################
    # print "NDleafs: ", len(NDLeafs)
    # Keep track of evals
    numEval = 0
    for leaf in NDLeafs:
        #print leaf.east.evaluated, leaf.east.cells
        #print leaf.west.evaluated, leaf.west.cells
        #exit()
        #Pull out the Poles
        east = leaf.table.rows[0]
        west = leaf.table.rows[-1]
        #Evaluate those poles if needed
        if not east.evaluated:
            for o,objScore in enumerate(problem.evaluate(east.cells)):
                east.cells[-(len(problem.objectives)-o)] = objScore
            east.evaluated = True
            numEval += 1
        if not west.evaluated:
            for o,objScore in enumerate(problem.evaluate(west.cells)):
                west.cells[-(len(problem.objectives)-o)] = objScore
            west.evaluated = True
            numEval += 1

        #Score the poles
        n = len(problem.decisions)
        weights = []
        for obj in problem.objectives:
              # w is negative when we are maximizing that objective
              if obj.lismore:
                  weights.append(+1)
              else:
                  weights.append(-1)
        weightedWest = [c*w for c,w in zip(west.cells[n:], weights)]
        weightedEast = [c*w for c,w in zip(east.cells[n:], weights)]
        westLoss = loss(weightedWest, weightedEast, mins = [obj.low for obj in problem.objectives], maxs = [obj.up for obj in problem.objectives])
        eastLoss = loss(weightedEast, weightedWest, mins = [obj.low for obj in problem.objectives], maxs = [obj.up for obj in problem.objectives])

        #Determine better Pole
        if eastLoss < westLoss:
            SouthPole, NorthPole = east, west
        else:
            SouthPole, NorthPole = west, east

        #Magnitude of the mutations
        g = abs(SouthPole.x - NorthPole.x)

        # print "leaf table rows: ", len(leaf.table.rows)


        #Iterate over the individuals of the leaf
        for row in leaf.table.rows:
            # print "row cells: ", len(row.cells)

            #Make a copy of the row in case we reject it
            copy = [item for item in row.cells]
            cx   = row.x

            for attr in range(0, len(problem.decisions)):

                #just some naming shortcuts
                me   = row.cells[attr]
                good = SouthPole.cells[attr]
                bad  = NorthPole.cells[attr]
                dec  = problem.decisions[attr]
                #print "dec: ", dec


                #Find direction to mutate (Want to mutate towards good pole)
                if me > good:  d = -1
                if me < good:  d = +1
                if me == good: d =  0

                row.cells[attr] = min(dec.up, max(dec.low, me + me*g*d))

            #Project the Mutant
            a    = row.distance(NorthPole)
            b    = row.distance(SouthPole)
            c    = NorthPole.distance(SouthPole)
            x    = (a**2 + row.c**2 - b**2) / (2*row.c+0.00001)  # this should not be row.c but c

            #Test Mutant for Acceptance
            GAMMA = 0.15 #note: make this a property #Vivek: I think this should not be here

            #print abs(cx-x), (cx + (g * GAMMA))
            if abs(x-cx) > (g * GAMMA) or problem.evalConstraints(row.cells[:n]): #reject it
                row.cells = copy
                row.x = x


    # After mutation; Convert back to JMOO Data Structures
    population = []
    for leaf in NDLeafs:
        for row in leaf.table.rows:
            if row.evaluated:
                population.append(jmoo_individual(problem, [x for x in row.cells[:len(problem.decisions)]], [x for x in row.cells[len(problem.decisions):]]))
            else:
                population.append(jmoo_individual(problem, [x for x in row.cells[:len(problem.decisions)]], None))

    # Return selectees and number of evaluations
    return population, numEval

def galeMutate2(problem, NDLeafs):

    #################
    # Mutation Phase
    #################
    print "NDleafs: ", len(NDLeafs)
    # Keep track of evals
    numEval = 0
    for leaf in NDLeafs:
        #print leaf.east.evaluated, leaf.east.cells
        #print leaf.west.evaluated, leaf.west.cells
        #exit()
        #Pull out the Poles
        east = leaf.table.rows[0]
        west = leaf.table.rows[-1]
        #Evaluate those poles if needed
        if not east.evaluated:
            for o,objScore in enumerate(problem.evaluate(east.cells)):
                east.cells[-(len(problem.objectives)-o)] = objScore
            east.evaluated = True
            numEval += 1
        if not west.evaluated:
            for o,objScore in enumerate(problem.evaluate(west.cells)):
                west.cells[-(len(problem.objectives)-o)] = objScore
            west.evaluated = True
            numEval += 1

        if not leaf.north.evaluated:
            for o,objScore in enumerate(problem.evaluate(leaf.north.cells)):
                leaf.north.cells[-(len(problem.objectives)-o)] = objScore
            leaf.north.evaluated = True
            numEval += 1

        if not leaf.south.evaluated:
            for o,objScore in enumerate(problem.evaluate(leaf.south.cells)):
                leaf.south.cells[-(len(problem.objectives)-o)] = objScore
            leaf.south.evaluated = True
            numEval += 1


        #Score the poles
        n = len(problem.decisions)
        weights = []
        for obj in problem.objectives:
              # w is negative when we are maximizing that objective
              if obj.lismore:
                  weights.append(+1)
              else:
                  weights.append(-1)
        weightedWest = [c*w for c,w in zip(west.cells[n:], weights)]
        weightedEast = [c*w for c,w in zip(east.cells[n:], weights)]
        weightedNorth = [c*w for c,w in zip(leaf.north.cells[n:], weights)]
        weightedSouth = [c*w for c,w in zip(leaf.south.cells[n:], weights)]


        # print leaf.south.id, leaf.north.id, leaf.east.id, leaf.west.id

        westLoss = loss(weightedWest, weightedEast, mins = [obj.low for obj in problem.objectives], maxs = [obj.up for obj in problem.objectives])
        eastLoss = loss(weightedEast, weightedWest, mins = [obj.low for obj in problem.objectives], maxs = [obj.up for obj in problem.objectives])
        northLoss = loss(weightedNorth, weightedSouth, mins = [obj.low for obj in problem.objectives], maxs = [obj.up for obj in problem.objectives])
        southLoss = loss(weightedSouth, weightedNorth, mins = [obj.low for obj in problem.objectives], maxs = [obj.up for obj in problem.objectives])

        westeast = abs(westLoss - eastLoss)
        northsouth = abs(northLoss - southLoss)

        # print westeast, northsouth


        if westeast > northsouth: # improvement is more in the west east direction
            #Determine better Pole
            if eastLoss < westLoss:
                SouthPole, NorthPole = east, west
            else:
                SouthPole, NorthPole = west, east
        else:

            if northLoss < southLoss:
                SouthPole, NorthPole = leaf.north, leaf.south
            else:
                SouthPole, NorthPole = leaf.south, leaf.north


        #Magnitude of the mutations
        g = abs(SouthPole.x - NorthPole.x)

        # print "leaf table rows: ", len(leaf.table.rows)


        #Iterate over the individuals of the leaf
        for row in leaf.table.rows:
            # print "row cells: ", len(row.cells)

            #Make a copy of the row in case we reject it
            copy = [item for item in row.cells]
            cx   = row.x

            for attr in range(0, len(problem.decisions)):

                #just some naming shortcuts
                me   = row.cells[attr]
                good = SouthPole.cells[attr]
                bad  = NorthPole.cells[attr]
                dec  = problem.decisions[attr]
                #print "dec: ", dec


                #Find direction to mutate (Want to mutate towards good pole)
                if me > good:  d = -1
                if me < good:  d = +1
                if me == good: d =  0

                row.cells[attr] = min(dec.up, max(dec.low, me + me*g*d))

            #Project the Mutant
            a    = row.distance(NorthPole)
            b    = row.distance(SouthPole)
            c    = NorthPole.distance(SouthPole)
            x    = (a**2 + c**2 - b**2) / (2*c+0.00001)
            y =  abs(a**2-(x/(2*c))**2)**0.5

            #Test Mutant for Acceptance
            GAMMA = 0.15 #note: make this a property #Vivek: I think this should not be here

            #print abs(cx-x), (cx + (g * GAMMA))
            if abs(x-cx) > (g * GAMMA) or problem.evalConstraints(row.cells[:n]): #reject it
                row.cells = copy
                row.x = x


    # After mutation; Convert back to JMOO Data Structures
    population = []
    for leaf in NDLeafs:
        for row in leaf.table.rows:
            if row.evaluated:
                population.append(jmoo_individual(problem, [x for x in row.cells[:len(problem.decisions)]], [x for x in row.cells[len(problem.decisions):]]))
            else:
                population.append(jmoo_individual(problem, [x for x in row.cells[:len(problem.decisions)]], None))

    # Return selectees and number of evaluations
    return population, numEval

def galeRegen(problem, unusedSlot, mutants, MU):
    
    howMany = MU - len(mutants)
    #print "HowMany: ",howMany
    # exit()
    
    # Generate random individuals
    population = []
    for i in range(howMany):
        population.append(jmoo_individual(problem, problem.generateInput(), None))
    
    return mutants+population, 0


def galeRegen2(problem, unusedSlot, mutants, MU):
    def trim(mutated, low, up):
        return max(low, min(mutated, up))

    def three_others(individuals, one):
        seen = [one]

        def other():
            while True:
                random_selection = random.randint(0, len(individuals) - 1)
                if individuals[random_selection] not in seen:
                    seen.append(individuals[random_selection])
                    break
            return individuals[random_selection]

        return other(), other(), other()

    def extrapolate(problem, individuals, f, cf):
        # #print "Extrapolate"
        one = individuals[random.randint(0, len(individuals)-1)]
        two, three, four = three_others(individuals, one)
        solution = []
        for d, decision in enumerate(problem.decisions):
            assert isinstance(two, jmoo_individual)
            x, y, z = two.decisionValues[d], three.decisionValues[d], four.decisionValues[d]
            if random.random() < cf:
                mutated = x + f * (y - z)
                solution.append(trim(mutated, decision.low, decision.up))
            else:
                solution.append(one.decisionValues[d])
        return jmoo_individual(problem, [float(d) for d in solution], None)


    howMany = MU - len(mutants)

    # Generate random individuals
    population = []
    for i in range(howMany):
        population.append(extrapolate(problem, mutants, jmoo_properties.F, jmoo_properties.CF ))

    return mutants+population, 0