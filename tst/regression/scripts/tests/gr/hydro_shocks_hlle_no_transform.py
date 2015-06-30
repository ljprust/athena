# Test script for relativistic hydro shock tubes in GR with HLLE without transforming

# Modules
import numpy as np
import scripts.utils.athena as athena
import scripts.utils.comparison as comparison

# Prepare Athena++
def prepare():
  athena.configure('g',
      prob='shock_tube_rel',
      coord='minkowski')
  athena.make()

# Run Athena++
def run():
  arguments = [
      '',
      'output1/file_type=vtk',
      'output1/variable=cons',
      'output1/dt=0.4',
      'time/cfl_number=0.4',
      'time/tlim=0.4',
      'mesh/nx1=400']
  for i in range(1,5):
    arguments[0] = 'job/problem_id=gr_hydro_shock' + repr(i)
    athena.run('hydro_sr/athinput.mb_'+repr(i), arguments)

# Analyze outputs
def analyze():
  headers = [('dens',), ('Etot',), ('mom',0)]
  tols = [[0.02,0.01,0.01], [0.01,0.01,0.02], [0.01,0.01,0.02], [0.5,0.01,0.02]]
  for i in range(1,5):
    x_ref,_,_,data_ref = athena.read_vtk('data/sr_hydro_shock{0}_hllc.vtk'.format(i))
    x_new,_,_,data_new = \
        athena.read_vtk('bin/gr_hydro_shock{0}.block0.out1.00001.vtk'.format(i))
    tols_particular = tols[i-1]
    for header,tol in zip(headers,tols_particular):
      array_ref = data_ref[header[0]]
      array_ref = array_ref[0,0,:] if len(header) == 1 else array_ref[0,0,:,header[1]]
      array_new = data_new[header[0]]
      array_new = array_new[0,0,:] if len(header) == 1 else array_new[0,0,:,header[1]]
      if header[0] == 'Etot':
        array_new = -array_new   # sign difference between SR and GR
      eps = comparison.l1_diff(x_ref, array_ref, x_new, array_new)
      eps /= comparison.l1_norm(x_ref, array_ref)
      if eps > tol or np.isnan(eps):
        return False
  return True