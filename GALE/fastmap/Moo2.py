

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

from BinaryTree import *
from The import *
from lib import *
from Row import *
from utility import *
import random
import math


class Moo2(BinaryTree):
    #t: Table()
    #big_n: len(t.rows)

    def __init__(i, problem, t, big_n, N=1):  # problem, table, number of rows in the table, N?

        BinaryTree.__init__(i)
        i.big_n = big_n
        i.table = t
        i.abort=False
        i.abortRate = 0.0
        i.east, i.west, i.c, i.x, i.north, i.south = None,None,None,None,None,None
        i.N = N
        i.problem = problem

    def project(i,rows):
        "Uses the O(2N) Fastmap heuristic."
        w    = one(rows) # any row, selected at random
        i.west = w.furthest()
        i.east = i.west.furthest()
        i.c    = i.west.distance(i.east)
        for row in rows:
            a    = row.distance(i.west)
            b    = row.distance(i.east)
            row.x= (a**2 + i.c**2 - b**2) / (2*i.c+0.00001)
            row.y= abs(a**2-(row.x/(2*i.c))**2)**0.5
            row.c = i.c
            row.a = a
            row.b = b
        rows = sorted(rows, key= lambda row: row.x)
        return rows

    def project2(i, rows):
        "Uses the O(2N) Fastmap heuristic."
        rows = sorted(rows, key= lambda row: row.y)
        #print len(rows)
        #for row in rows: print row.id, row.x, row.y
        i.north = rows[0]
        i.south = rows[-1]
        return rows

    def split(i,rows,mid,parent):
        i.x = rows[mid].x
        i.parent = parent
        i.east = rows[0]
        i.west = rows[-1]
        return rows[:mid],rows[mid:]

    def split2(i,rows,mid,parent):
        i.y = rows[mid].y
        i.parent = parent
        i.north = rows[0]
        i.south = rows[-1]
        return rows[:mid], rows[mid:]

    def divide(i,abort=False,minnie=30):  # minnie: this is stopping criteria (N ** 0.5)
        def aFew(rows) :
            all = map(lambda r: r.cells,rows)
            # print "# of rows: ", len(rows), rows[0]
            # print "# of alls: ", len(all), all[0]
            # exit()
            new_t = i.table.clone(some(all,The.alpha))
            for ind,irow in enumerate(new_t.rows):
                irow.evaluated = rows[ind].evaluated
            return new_t

            # Put objectives into its own array
        pop = []
        # print "Length of the table is : ", len(i.table.rows)  #Table -> row -> cell

        for row in i.table.rows:
            # print "Objectives of table is : ", row.table.objectives
            # row.table is like a backward loop to get objectives etc
            pop.append([x for obj, x, y in all(row.table.objectives, row, row)])



        # grab number of rows
        n = len(i.table.rows)

        # Project the rows onto 1D
        # now i.table.rows is sorted according to the x from fast map
        i.table.rows = i.project(i.table.rows)
        # It is correct
        # print len(i.table.rows)
        # import time
        # time.sleep(1)
        test = i.project2(i.table.rows)


        # Find a good splitting point
        m, _ = i.binaryChop(i.table.rows, n/2, None, 2*n ** 0.5, n)

        #Vivek: m is the cut which achieves min(lhsdelta - rhsdelta)

        # Proceed if not aborted
        i.abort = The.allowDomination and abort
        if not i.abort and n >= minnie :

            # Do the Split
            wests, easts = i.split(i.table.rows, m, n)
            i.split2(i.table.rows, 0, 0)


            # Precautions if too many reps in east & west
            if i.west != i.east:
                if (i.N > m):
                    littleN = m
                else:
                    littleN = i.N

                #print "littleN: ", littleN

                westAbort = False
                eastAbort = False
                if not i.east.evaluated:
                    i.east.evaluated = True
                    for o,objScore in enumerate(i.problem.evaluate(i.east.cells)):
                        i.east.cells[-(len(i.problem.objectives)-o)] = objScore
                if not i.west.evaluated:
                    i.west.evaluated = True
                    for o,objScore in enumerate(i.problem.evaluate(i.west.cells)):
                        i.west.cells[-(len(i.problem.objectives)-o)] = objScore

                # print "East: ", i.east
                weights = []
                for obj in i.problem.objectives:
                    # w is negative when we are maximizing that objective
                    if obj.lismore:
                        weights.append(+1)
                    else:
                        weights.append(-1)
                k = len(i.problem.objectives)
                weightedWest = [c*w for c,w in zip(i.west.cells[-k:], weights)]
                weightedEast = [c*w for c,w in zip(i.east.cells[-k:], weights)]
                westLoss = loss(weightedWest, weightedEast, mins = [obj.low for obj in i.problem.objectives], maxs = [obj.up for obj in i.problem.objectives])
                eastLoss = loss(weightedEast, weightedWest, mins = [obj.low for obj in i.problem.objectives], maxs = [obj.up for obj in i.problem.objectives])
                EPSILON = 1.0
                # print "East: ", i.east
                # print "West: ", i.west

                if westLoss < EPSILON * eastLoss:
                    eastAbort = True    # east dominates
                if eastLoss < EPSILON * westLoss:
                    westAbort = True   # west dominates
                # print eastAbort, westAbort

                # Copy into each "half"
                i.lhs = Moo2(i.problem, aFew(wests), i.big_n, littleN)
                i.rhs = Moo2(i.problem, aFew(easts), i.big_n, littleN)

                # Divide each "half"
                i.rhs.divide(abort=eastAbort, minnie=minnie)
                i.lhs.divide(abort=westAbort, minnie=minnie)


        return i

    def binaryChop(i, rows, cut, delta, min_n, lastcut=None):
        "perform binary chop to find an appropriate place to split"
        # print "min_n: ", min_n

        # stop if too small
        if cut < min_n or lastcut-cut < min_n:
            return cut,delta

        # segment left and right sides
        left = rows[:cut]
        right = rows[cut:]

        # print "Length of left: ", len(left),
        # print " Length of right: ", len(right)

        # get spreads (VARIANCE) of each side
        z = len(i.problem.decisions)
        leftSpread = spacing([l.cells[:z] for l in left])  # var([l.x for l in left])
        rightSpread = spacing([r.cells[:z] for r in right])  # var([r.x for r in right])
        delta = abs(leftSpread - rightSpread)  # sdiv
        # print "Cut: ", cut
        # print delta, cut
        #recurse
        lhscut,lhsdelta = i.binaryChop(rows, cut/2, delta, min_n, cut)
        rhscut,rhsdelta = i.binaryChop(rows, cut + (lastcut - cut)/2, delta, min_n, cut)

        #minimize deltas
        smallest = min(delta, lhsdelta, rhsdelta)
        if smallest == delta:
            # print "Orginal: ", cut
            return cut, delta
        elif smallest == lhsdelta:
            # print "LHS: ", lhscut
            return lhscut, lhsdelta
        else:
            # print "RHS: ", rhscut
            return rhscut, rhsdelta

    def __repr__(i):
        n = "%1.0f" % (100.0*((len(i.table.rows)) / (float)(i.big_n)))
        n += "%"
        s = i.table.status()
        post = "<-- pruned (" + str(i.abortRate) + ")" if i.abort else "(" + str(i.abortRate) + ")"
        if i.c:
            return "#%s / %4.4s = %4.4s : %s" % (n,i.c,i.x,post)
        else:
            return "#%s %s %s" % (n,s,post)
