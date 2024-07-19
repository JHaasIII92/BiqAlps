#! /usr/bin/env python
################################################################################
import numpy as np
import gurobipy as grb
import sys

################################################################################

if len(sys.argv) != 2:
    print('Usage:  ./biq.py FILENAME')

else:
    # Get filename
    FILENAME = sys.argv[1]
    print('Using file:', FILENAME)

    # Load the graph
    data = np.loadtxt(FILENAME, skiprows=1, dtype=int)
    n = data[:,:2].max()

    # Initialize Model
    m = grb.Model('biq')

    # Introduce Variables
    print('Introducing variables...')
    x = {}
    for i in range(n):
        x[i] = m.addVar(vtype=grb.GRB.BINARY, name='x%d' % (i+1))
    m.update()

    # Define Objective
    print('Defining the objective function...')
    obj = 0
    for i, j, val in data:
        obj += val*x[i-1]*x[j-1]
        if i != j:
            obj += val*x[j-1]*x[i-1]
    m.setObjective(obj, sense=grb.GRB.MINIMIZE)

    # Optimize the Model
    print('============')
    m.setParam('TimeLimit', 600)
    m.optimize()
    print('============')

    # print the objective value
    print('Minimum objective = {0:.0f}'.format(obj.getValue()))

    # print the solution
    outstring = 'X = {'
    for i in range(n):
        if x[i].x == 1.0:
            outstring += '{0:>3} '.format(i+1)
    outstring += '}'
    print(outstring)
