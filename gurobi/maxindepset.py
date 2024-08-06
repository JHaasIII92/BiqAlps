#! /usr/bin/env python3
################################################################################
import numpy as np
import gurobipy as grb
import sys

################################################################################

if len(sys.argv) != 2:
    print('Usage:  ./maxindepset.py FILENAME')

else:
    # Get filename
    FILENAME = sys.argv[1]
    print('Using file:', FILENAME)

    # Load the graph
    data = np.loadtxt(FILENAME)
    n = int(data[0,0])

    # Initialize Model
    m = grb.Model('max-indep-set')

    # Introduce Variables
    print('Introducing variables...')
    x = {}
    for i in range(n):
        x[i] = m.addVar(vtype=grb.GRB.BINARY, name='x%d' % (i+1))
    m.update()

    # Define Objective
    print('Defining the objective function...')
    obj = 0
    for i in range(n):
        obj += x[i]

    m.setObjective(obj, sense=grb.GRB.MAXIMIZE)

    # Define Constraints
    k = 0
    for i in range(1, n+1):
        for j in range(1, i):
            foundij = 0
            for ii, jj in data[1:,:]:
                if i == ii and j == jj:
                    foundij = 1
            if foundij == 0:
                m.addQConstr(x[i-1]*x[j-1], grb.GRB.EQUAL, 0, 'c%d' % k)
                k = k + 1
                print(k, i, j)

    # Optimize the Model
    print('============')
    m.setParam('TimeLimit', 600)
    m.optimize()
    print('============')

    # Print the objective value
    print('Max-Independent-Set = {0:.0f}'.format(obj.getValue()))

    # Print the solution
    outstring = 'X = {'
    for i in range(n):
        if x[i].x == 1.0:
            outstring += '{0:>3} '.format(i+1)
    outstring += '}'
    print(outstring)
