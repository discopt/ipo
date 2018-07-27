#!/usr/bin/python

import sys
import IPO
import random

#objectiveType = 'weights'
#objectiveType = 'random'
#objectiveType = '0/1-random'
objectiveType = '0/-1-random'
objectiveZeroProbability = 0.5

N = 10000

V = range(10)
E = [ (i,j) for i in xrange(len(V)) for j in xrange(i+1,len(V)) ]

# Write LP.

fileName = '/tmp/ipo.lp'
lpData = '''
minimize cost: -x#0#1
s.t.
'''
for v in V:
  lpData += '  ' + ' + '.join([ 'x#%d#%d' % (i,j) for (i,j) in E if i == v or j == v ]) + ' == 1\n'
lpData += 'generals\n'
lpData += 'binaries\n'
for e in E:
  lpData += '  x#%d#%d\n' % e
lpData += 'end\n'
lpFile = open(fileName, 'w')
lpFile.write(lpData)
lpFile.close()

mils = IPO.MixedIntegerLinearSet(fileName)
matchingOracle = IPO.SCIPOracle(fileName, mils)
dominantOracle = IPO.DominantOracle("Dom", matchingOracle)
oracle = dominantOracle
P = IPO.Polyhedron(oracle)

cons = mils.getConstraints()
lp = IPO.LinearProgram(P.space)

dim, (points, rays), eqns = IPO.affineHull(oracle)

print 'Generating facets for %d directions.' % (N)
for i in xrange(N):
#  print('Random objective #%d' % (i+1))
  objective = [None] * lp.numVariables
  if objectiveType == 'weights':
    for c in xrange(lp.numVariables):
      objective[c] = weights[c] + random.random()
  elif objectiveType == 'random':
    for c in xrange(lp.numVariables):
      objective[c] = random.random() - 0.5
  elif objectiveType == '0/1-random':
    for c in xrange(lp.numVariables):
      rnd = random.random()
      if rnd <= objectiveZeroProbability:
        objective[c] = 0
      else:
        objective[c] = 1
  elif objectiveType == '0/-1-random':
    for c in xrange(lp.numVariables):
      rnd = random.random()
      if rnd <= objectiveZeroProbability:
        objective[c] = 0
      else:
        objective[c] = -1
  elif objectiveType == '-1/0/1-random':
    for c in xrange(lp.numVariables):
      if random.random() <= objectiveZeroProbability:
        objective[c] = 0
      elif random.random() > 0.5:
        objective[c] = 1
      else:
        objective[c] = -1
  else:
    assert False
  lp.changeObjective(objective)

  while True:
#    sys.stdout.write('.')
#    sys.stdout.flush()
#    print('Solving LP with %d rows.' % (lp.numRows))
    (status, vector, value) = lp.solve()
    if status == IPO.LinearProgram.OPTIMAL:
#      print('LP-optimum: %s' % (vector))
      con = P.separatePoint(vector)
      if con is None:
        break
      print('Separated facet: %s' % (con))
      lp.addConstraint(con)
    elif status == IPO.LinearProgram.UNBOUNDED:
#      print('LP-ray: %s' % (vector))
      con = P.separateRay(vector)
      if con is None:
        break
      print('Separated facet: %s' % (con))
      lp.addConstraint(con)
    else:
      raise Exception('Unexpected status %d' % (status))
print
