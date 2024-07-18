#! /usr/bin/env python3
################################################################################
import numpy as np
import gurobipy as grb
import sys

################################################################################

if len(sys.argv) != 2:
    print('Usage:  ./kcluster.py FILENAME')

else:
    # Get filename
    FILENAME = sys.argv[1]
    print('Using file:', FILENAME)

    # Load the graph
    data = np.loadtxt(FILENAME)
    n = int(data[0,0])
    k = int(data[0,1])
    #n = data[:,:2].max()

    # Initialize Model
    m = grb.Model('kcluster')

    # Introduce Variables
    print('Introducing variables...')
    x = {}
    for i in range(n):
        x[i] = m.addVar(vtype=grb.GRB.BINARY, name='x%d' % (i+1))
    m.update()

    # Define Objective
    print('Defining the objective function...')
    obj = 0
    for i, j, val in data[1:,:]:
        obj += val*x[i-1]*x[j-1]
        obj += val*x[j-1]*x[i-1]

    m.setObjective(obj, sense=grb.GRB.MAXIMIZE)

    # Define Constraint
    sumx = 0
    for i in range(n):
        sumx += x[i]
    m.addLConstr(sumx, grb.GRB.EQUAL, k, "c0")

    # Optimize the Model
    print('============')
    m.setParam('TimeLimit', 600)
    m.optimize()
    print('============')

    # print the objective value
    print('Max-k-cluster = {0:.0f}'.format(obj.getValue()))

    # print the solution
    outstring = 'X = {'
    for i in range(n):
        if x[i].x == 1.0:
            outstring += '{0:>3} '.format(i+1)
    outstring += '}'
    print(outstring)
