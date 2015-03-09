import os
import sys
import inspect
import random

from jmoo_individual import *
from jmoo_algorithms import *
from jmoo_stats_box import *


currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir)

import jmoo_properties


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


def trim(mutated, low, up):
    return max(low, min(mutated, up))


def extrapolate(problem, individuals, one, f, cf):
    # #print "Extrapolate"
    two, three, four = three_others(individuals, one)
    # #print two,three,four
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

def better(problem,individual,mutant):
    # From Joe: Score the poles
    n = len(problem.decisions)
    weights = []
    for obj in problem.objectives:
        # w is negative when we are maximizing that objective
        if obj.lismore:
            weights.append(+1)
        else:
            weights.append(-1)
    weighted_individual = [c*w for c,w in zip(individual.fitness.fitness, weights)]
    weighted_mutant = [c*w for c,w in zip(mutant.fitness.fitness, weights)]
    individual_loss = loss(weighted_individual, weighted_mutant, mins = [obj.low for obj in problem.objectives], maxs = [obj.up for obj in problem.objectives])
    mutant_loss = loss(weighted_mutant, weighted_individual, mins = [obj.low for obj in problem.objectives], maxs = [obj.up for obj in problem.objectives])

    if individual_loss < mutant_loss:
        return mutant
    else:
        return individual  # otherwise


def de_selector(problem, individuals):
    #print "selector"
    newer_generation = []
    for individual in individuals:
        if not individual.valid:
            individual.evaluate()
    no_evals = len(individuals)
    #print "Length of population: ", len(individuals)
    #print "F: ", jmoo_properties.F
    #print "CF: ", jmoo_properties.CF
    for individual in individuals:
        #print "Old Decision: ",individual.decisionValues
        #print "Old Score: ", individual.fitness.fitness
        mutant = extrapolate(problem, individuals, individual, jmoo_properties.F, jmoo_properties.CF)
        mutant.evaluate()
        no_evals += 1
        ##print "Mutant Decisions: ",mutant.decisionValues
        #print "New Score: ", mutant.fitness.fitness
        newer_generation.append(better(problem, individual, mutant))

    #print len(newer_generation)
    return newer_generation, no_evals

#Vivek: This is just a stub
def de_mutate(problem, population):
    #print "mutate"
    #print "Length of the Population: ",len(population)
    #for i,p in enumerate(population):
    #    print ">>", i, sum(p.fitness.fitness)
    #print "Minimum in population: ", min([sum(p.fitness.fitness) for p in population])

    return population, 0

#Vivek: This is just a stub
def de_recombine(problem, unusedSlot, mutants, MU):
    #print "recombine"
    return mutants, 0