import random
import IPO
import sys

if len(sys.argv) > 2:
  numRandom = int(sys.argv[2])
else:
  numRandom = 0

(P,mils,origObj) = IPO.readSCIP(sys.argv[1])

print('Dimension: %d' % (P.dimension))

for i in xrange(-1, numRandom):
  if i < 0:
    print('\nOriginal objective')
    objective = origObj
  else:
    print('\nRandom objective #%d' % (i))
    objective = [ (random.random() - 0.5) for j in xrange(P.space.dimension) ]
  for facet in IPO.generateFacets(P, mils, objective):
    print facet
    mils.addConstraint(facet)
