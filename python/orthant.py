#!/usr/bin/python -O

import ipo

def init(tempDirectory, instance):
  try:
    n = int(instance[0])
  except:
    raise Exception("Invalid instance \"%s\"." % (' '.join(instance)))
  return { 'variables' : [ 'x#%d' % (v,) for v in xrange(n) ] }

def optimize(tempDirectory, instance, objective, forceOptimal):
  try:
    n = int(instance[0])
  except:
    raise Exception("Invalid instance \"%s\"." % (' '.join(instance)))
  solution = [ 0 for v in xrange(n) ]
  for v in xrange(n):
    if objective[v] > 0:
      solution[v] = 1 # Create v'th unit vector.
      return ('unbounded', solution)
  return ('optimal', solution)

ipo.oracleCall(init, optimize)
