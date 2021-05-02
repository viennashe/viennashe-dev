#!/usr/bin/env python
##============================================================================
##   Copyright (c) 2011-2014, Institute for Microelectronics,
##                            Institute for Analysis and Scientific Computing,
##                            TU Wien.
##
##                            -----------------
##     ViennaSHE - The Vienna Spherical Harmonics Expansion Boltzmann Solver
##                            -----------------
##
##                    http://viennashe.sourceforge.net/
##
##   License:         MIT (X11), see file LICENSE in the base directory
##===============================================================================

import sys
import os


## @example resistor.py A simple 1D resistor simulation using the Python interface

# append the path of the viennashe python bindings
sys.path.append('../../pythonbinding/src')
sys.path.append('./pythonbinding/src')

try:
  import pyviennashe as viennashe
except Exception, ex:
  print ex
  print "\nCannot import pyviennashe. Possible causes:"
  print "\n - The LD_LIBRARY_PATH (UNIX) or the WINDOWS paths are set incorrectly."
  print "\n - The PYTHONPATH needs to be set to the same directory as the library search path (LD_LIBRARY_PATH for *NIX SYSTEMS)"
  sys.exit(1)

viennashe.initalize()


len_x = 1e-6
points_x = 51

SI_ID    = viennashe.get_silicon_id()[1]
METAL_ID = viennashe.get_metal_id()[1]
SIO2_ID  = viennashe.get_sio2_id()[1]
HFO2_ID  = viennashe.get_hfo2_id()[1]

dev = viennashe.create_1d_device(len_x, points_x)

num_cells = viennashe.get_num_cells(dev)[1]

print "num_vertices = ", viennashe.get_num_vertices(dev)[1]
print "num_cells    = ", num_cells


matids=[]
for i in range(0, num_cells): matids.append(SI_ID)
matids[0] = METAL_ID
matids[num_cells-1] = METAL_ID

Nd = []; Na = []
for i in range(0, num_cells):
  Nd.append(1e25)
  Na.append(1e7)

bnd_cells = [0, num_cells-1]
bnd_pot = [0.0, 1.0]

viennashe.initalize_device(dev, matids, Nd, Na)
viennashe.set_contact_potential_cells(dev, bnd_cells, bnd_pot, 2)


## Create the simulator configuration
conf_dd = viennashe.create_config();

## Apply standard DD configuration
viennashe.config_standard_dd(conf_dd)
viennashe.set_nonlinear_solver_config(conf_dd, viennashe.nonlinear_solver_gummel, 55, 0.5)

## Get the simulator instance
sim_dd = viennashe.create_simulator(dev, conf_dd)

## This is DD, so we do not need to give any inital guess
viennashe.run(sim_dd)

reg = viennashe.create_quantity_register(sim_dd )

viennashe.write_to_gnuplot(reg, "Electrostatic potential",     "potential_she.dat")
viennashe.write_to_gnuplot(reg, "Electron density",            "n_she.dat")
viennashe.write_to_gnuplot(reg, "Hole density",                "p_she.dat")

## FIN
viennashe.finalize()

