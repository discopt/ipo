import sys
import os
import gurobipy
from gurobipy import quicksum, GRB
import random
from xml.dom import minidom
import gzip

MAX_NUM_VARS = 12000
TIME_LIMIT = 600
ATTEMPTS = 10
MIN_SUCCESSFUL = 3

env = gurobipy.Env(params={'outputFlag':0})

def sampleObjective(dimension):
  squaredNorm = 0.0
  while squaredNorm < 1.0e-5:
    vector = [ random.gauss(0.0, 1.0) for _ in range(dimension) ]
    squaredNorm = sum( x**2 for x in vector )
  norm = squaredNorm**0.5
  return tuple([ x/norm for x in vector ])


for instance_file in sys.argv[1:]:
  instance_name = instance_file.rsplit(os.sep, 1)[1].split('.', 1)[0]

  model = gurobipy.read(instance_file, env)
  if model.numVars > MAX_NUM_VARS:
    continue

  variables = model.getVars()
  scale = 1 if model.modelSense == GRB.MAXIMIZE else -1
  model_objective = { var.varName: scale * var.obj for var in variables }
  model.modelSense = GRB.MAXIMIZE
  model.params.timeLimit = TIME_LIMIT

  print(f'Instance "{instance_name}" has {model.numVars} variables.')
  generating_objectives = []
  target_solutions = []
  for attempt in range(1, ATTEMPTS+1):
    random_objective = sampleObjective(len(variables))
#    for i,var in enumerate(variables):
#      print(f'  {var} -> {model_objective[i]} {random_objective[i]}')

    model.setObjective(gurobipy.quicksum( random_objective[i] * var for i,var in enumerate(variables) ) )
    model.optimize()

    print(f'  Attempt #{attempt} ran {model.runtime:.1f}s with status {model.status}.')

    if model.status == GRB.OPTIMAL:
      generating_objectives.append( { var.varName: random_objective[i] for i,var in enumerate(variables) } )
      solution = { var.varName: var.x for var in variables }
      target_solutions.append( solution )

      if len(generating_objectives) >= MIN_SUCCESSFUL:
        break
    if model.status == GRB.INFEASIBLE:
      break

  if len(generating_objectives) >= MIN_SUCCESSFUL:
    for gen in range(len(generating_objectives)):
      target_objective = model_objective
      target_solution = target_solutions[gen]
      generating_objective = generating_objectives[gen]

      document = minidom.Document()
      xml_root = document.createElement('inverse-mip')
      xml_root.setAttribute('dimension', f'{len(variables)}')
      document.appendChild(xml_root)

      xml_target_objective = document.createElement('target-objective')
      xml_target_objective.setAttribute('norm', 'L1')
      for var in variables:
        xml_obj = document.createElement('variable')
        xml_obj.setAttribute('name', var.varName)
        xml_obj.setAttribute('coefficient', f'{target_objective[var.varName]}')
        xml_target_objective.appendChild(xml_obj)
      xml_root.appendChild(xml_target_objective)

      xml_mip = document.createElement('mip')
      xml_mip.setAttribute('sense', 'maximize')
      xml_mip.setAttribute('file', instance_file)
      xml_target_solution = document.createElement('target-solution')
      for var in variables:
        xml_sol = document.createElement('variable')
        xml_sol.setAttribute('name', var.varName)
        xml_sol.setAttribute('value', f'{target_solution[var.varName]}')
        xml_sol.setAttribute('gen-objective', f'{generating_objective[var.varName]}')
        xml_target_solution.appendChild(xml_sol)
      xml_mip.appendChild(xml_target_solution)
      xml_root.appendChild(xml_mip)

      gzip.open(f'{instance_name}_t{gen+1}.xml.gz', 'wt').write(xml_root.toprettyxml())

