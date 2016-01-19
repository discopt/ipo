#!/usr/bin/python -O

import ipo

# these modules are needed to call concorde
import math
import re
import os
import subprocess

# prefix for temporary files created to communicate with concorde
CONCORDE="/home/xammy/software/concorde/TSP/concorde"

def init(tempDirectory, instance):
  try:
    # read number of nodes
    nnodes = int(instance[0])
  except:
    raise Exception("Invalid instance \"%s\"." % (' '.join(instance)))
  # create oneavariable per undirected edge in the complete graph on 'nnodes' nodes
  return { 'variables' : [ 'x#%d#%d' % (i, j) for i in xrange(nnodes) for j in xrange(i + 1, nnodes) ] }

def optimize(tempDirectory, instance, objective, forceOptimal):
  try:
    nnodes = int(instance[0])
    nedges = nnodes * (nnodes - 1) / 2
  except:
    raise Exception("Invalid instance \"%s\"." % (' '.join(instance)))

  if nnodes <= 2:
    return ('infeasible', None)

  # Shift objective to be minimization of nonnegative vector.
  shift = max(objective)
  for i in xrange(nedges):
    objective[i] -= shift
  minObjective = min(objective)

  # Scale objective to be not too ugly.
  scale = 1
  while minObjective < -100000:
    scale *= 2
    minObjective /= 2
  for i in xrange(nedges):
    objective[i] /= scale

  # create .tsp-file to be read by concorde
  with open('%s/instance.tsp' % (tempDirectory), "w") as f:
    f.write("NAME : temporary\n")
    f.write("COMMENT : temporary tsp instance used for IPO oracle\n")
    f.write("TYPE : TSP\n")
    f.write("DIMENSION : %d\n" % nnodes)
    f.write("EDGE_WEIGHT_TYPE : EXPLICIT\n")
    f.write("EDGE_WEIGHT_FORMAT : UPPER_ROW\n")
    f.write("EDGE_WEIGHT_SECTION\n")
    for v in xrange(nedges):
      f.write("%d\n" % -objective[v])

  # run concorde
  with open('%s/instance.log' % (tempDirectory), "w") as logfile:
    p = subprocess.Popen([CONCORDE, 'instance.tsp'], stdout=logfile, stderr=logfile, cwd=tempDirectory)
    p.wait()
  if p.returncode:
    raise Exception("Error running concorde, see %s/instance.log." % tempDirectory)

  # read concorde's output as a permutation
  with open('%s/instance.sol' % (tempDirectory), "r") as f:
    permutation = [int(s) for s in f.read().split()[1:]]

  # translate permutation into a matrix
  solution_matrix = [nnodes * [0] for i in xrange(nnodes)]
  for i in xrange(nnodes):
    v = permutation[i]
    w = permutation[(i + 1) % nnodes]
    solution_matrix[v][w] = 1
    solution_matrix[w][v] = 1

  # construct solution vector from upper triangle of that matrix
  solution = []
  for i in xrange(nnodes):
    for j in xrange(i + 1, nnodes):
      solution.append(solution_matrix[i][j])

  # delete temporary files
  os.remove('%s/instance.tsp' % (tempDirectory))
  os.remove('%s/instance.log' % (tempDirectory))
  os.remove('%s/instance.sol' % (tempDirectory))

  return ('optimal', solution)

ipo.oracleCall(init, optimize)
