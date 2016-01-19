import sys
from string import split

def oracleCall(callbackInit, callbackMaximize):
  """
  Parses the arguments and calls one of the callbacks:
    
    callbackInit(instance): 
      This function is supposed to return a dict with these keys:
      'variables' followed by a list of d strings which are the variable names.

    callbackMaximize(instance, objective, forceOptimal):
      This function is supposed to return a pair of one of these types:
      ('infeasible', None)
      ('unbounded', direction) where direction is a list of d rationals representing an unbounded direction.
      ('optimal', solution) where solution is a list of d rationals representing an optimal point.
      ('feasible', solution) (only if forceOptimal=False) where solution is a list of d rationals representing a point.
  """

  if len(sys.argv) >= 3 and sys.argv[1] == '--init':
    for key,value in callbackInit(sys.argv[2], sys.argv[3:]).iteritems():
      if key == 'variables':
        print 'variables %d\n%s' % (len(value), '\n'.join(value))
      else:
        print 'ERROR %s' % (key,)
    return

  if len(sys.argv) >= 3 and sys.argv[1] == '--optimize':
    objective = map(int, split(sys.argv[2]))
    (status,vector) = callbackMaximize(sys.argv[3], sys.argv[4:], objective, True)
    if status == 'infeasible':
      print 'infeasible'
    elif status == 'unbounded' or status == 'optimal':
      print '%s\n%s' % (status, '\n'.join(map(str, vector)))
    else:
      print 'ERROR: %s' % (status,)
    return
  
  if len(sys.argv) >= 4 and sys.argv[1] == '--heuristic':
    objective = map(int, split(sys.argv[2]))
    (status,vector) = callbackMaximize(sys.argv[3], sys.argv[4:], objective, False)
    if status == 'infeasible':
      print 'infeasible'
    elif status == 'unbounded' or status == 'optimal' or status == 'feasible':
      print '%s\n%s' % (status, '\n'.join(map(str, vector)))
    else:
      print 'ERROR: %s' % (status,)
    return

  print 'An IPO oracle must be called in one of these ways:'
  print '%s --init INSTANCE...' % (sys.argv[0])
  print '%s --optimize OBJVEC INSTANCE...' % (sys.argv[0])
  print '%s --heuristic OBJVEC INSTANCE...' % (sys.argv[0])
  print
  print 'where OBJVEC is a single argument representing an integral vector with whitespace-separated coordinates,'
  print 'and where INSTANCE... specifies the instance using as many arguments as desired.'

