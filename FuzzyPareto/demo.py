import math, os, inspect, sys

currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir)

from jmoo_individual import *
from jmoo_properties import *


class LeftGaussianFunction:
    def __init__(i,individual):
        i.vector = individual
        i.c = -1
        i.sigma = 0.5 # This is heuristically chosen

    def do_it_x(i,x):
        if x <= i.c:
            return 1
        else:
            return math.exp(-0.5 * ((x-i.c)/i.sigma)**2)

    def do_it(i):
        return [i.do_it_x(vec) for vec in i.vector]


class DominanceRelation:
    def __init__(i, u, v):
        # u and v are two vectors of jmoo_individuals
        assert(u.fitness is not None), "u has not been evaluated"
        i.u = u.fitness.fitness
        assert(v.fitness is not None), "v has not been evaluated"
        i.v = v.fitness.fitness
        # u and v are only objectives, since domination only cares about objectives

        i.product_u = -1
        i.product_v = -1
        i.evaluate()

    def diff(i, u, v):
        """Finds the diff in the list"""
        assert(len(u) == len(v)), "The length of two vectors are not equal"
        return [a-b for a, b in zip(u, v)]

    def evaluate(i):
        # Performance
        perf_u = i.diff(i.u, i.v)
        perf_v = i.diff(i.v, i.u)

        # Gaussian Transform
        phi_u = LeftGaussianFunction(perf_u).do_it()
        phi_v = LeftGaussianFunction(perf_v).do_it()

        # Fuzzy Product Value: Product of all the elements of a list
        from operator import mul
        i.product_u = reduce(mul, phi_u, 1)
        i.product_v = reduce(mul, phi_v, 1)

    def domination_degree(i):
        temp_sum = i.product_u + i.product_v
        return i.product_u/temp_sum, i.product_v/temp_sum

def generate_population(problem,N):
    population = []
    for i in xrange(N):
        temp_decision = problem.generateInput()
        population.append(jmoo_individual(problem, temp_decision, problem.evaluate(temp_decision), index=i))
    return population

def fuzzy_fitness_measurement(individual, population):
    """
    :param individual: jmoo_individual
    :param population: list of jmoo_individual
    :return: FFM score of individual
    """
    FFM = 0
    for v in population:
        dom = DominanceRelation(individual,v)
        dom_individual, _ = dom.domination_degree()
        FFM += dom_individual

    return FFM/len(population)

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

    return False if individual_loss < mutant_loss else True

def fast_domination_sort(problem, population):
    non_dominated = []
    for p in population:
        count = 0
        for q in population:
            if p != q:  # don't compare p with q
                if better(problem, q, p):
                    count += 1
        if count == 0:
            non_dominated.append(p)
    return non_dominated



def baseline(problem, population):
    obj = len(population[0].fitness.fitness)
    up = [-1e6 for _ in xrange(obj)]
    low = [1e6 for _ in xrange(obj)]
    for p in population:
        for o in xrange(obj):
            if up[o] < p.fitness.fitness[o]: up[o] = p.fitness.fitness[o]
            elif low[o] > p.fitness.fitness[o]: low[o] = p.fitness.fitness[o]
    for o, obj in enumerate(problem.objectives):
        rangeX5 = up[o] - low[o]
        obj.low = -1*rangeX5
        obj.up = rangeX5
    return problem


def try_it():
    random.seed(3)
    problem = viennet2()
    population = generate_population(problem,100)
    problem = baseline(problem,population)
    for i in xrange(len(population)):
        individual = population[i]
        individual.FFM = fuzzy_fitness_measurement(individual, population[:i]+population[i:])

    for pop in fast_domination_sort(problem, population):
        print pop.FFM

    print max(pop.FFM for pop in population)

    what I am trying to do here is to find whether fuzzy pareto actuall works. The easiest way to check difference is
    to find a set of non dominated solutions and then find the fuzzy score.




if __name__ == "__main__":
    try_it()