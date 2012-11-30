#!/usr/bin/env python

import cubit
import boundary_definition
import cubit2specfem3d 

import os
import sys

xL, yL, zL = 3000.0, 3000.0, 1000.0 # zL refers to total length of SOLID, i.e. not including the water
water_thickness = 200.0
water_vp, water_vs, water_rho = 1480.0, 0.0, 1000.0

cubit.cmd('reset')
cubit.cmd('brick x '+str(xL)+' y '+str(yL)+' z '+str(water_thickness)) # centered at (0,0,0)
cubit.cmd('brick x '+str(xL)+' y '+str(yL)+' z '+str(zL)) # centered at (0,0,0)
center_xy = [xL/2, yL/2] # center point of domain; best to use a shot location here
cubit.cmd('volume 1 move x '+str(center_xy[0])+' y '+str(center_xy[1])+' z '+str(water_thickness/2.0)) # water layer
cubit.cmd('volume 2 move x '+str(center_xy[0])+' y '+str(center_xy[1])+' z -'+str(zL/2.0)) # solid
cubit.cmd('merge all')

# Meshing the volumes
elementsize = [100.0, 100.0]

cubit.cmd('volume 1 size '+str(elementsize[0]))
#cubit.cmd('mesh volume 1')
cubit.cmd('volume 2 size '+str(elementsize[1]))
cubit.cmd('mesh volume 1 2')

#### End of meshing 


node_list=cubit.parse_cubit_list('node','all')
print node_list[0]

###### This is boundary_definition.py of GEOCUBIT
#..... which extracts the bounding faces and defines them into blocks
boundary_definition.entities=['face']
boundary_definition.define_bc(boundary_definition.entities,parallel=True)

#### Define material properties for the 3 volumes ################
cubit.cmd('#### DEFINE MATERIAL PROPERTIES #######################')

cubit.cmd('block 1 name "acoustic 1" ')
cubit.cmd('block 1 attribute count 4')
cubit.cmd('block 1 attribute index 1 1')
cubit.cmd('block 1 attribute index 2 '+str(water_vp))   
cubit.cmd('block 1 attribute index 3 '+str(water_vs))   
cubit.cmd('block 1 attribute index 4 '+str(water_rho))  

cubit.cmd('block 2 name "acoustic tomography_model.xyz 1" ')        # elastic material region
cubit.cmd('block 2 attribute count 2')
cubit.cmd('block 2 attribute index 1 -1')      # flag for material: -1 for 1. undefined material
cubit.cmd('block 2 attribute index 2 2')      # flag for tomographic model 


#cubit.cmd('block 2 name "acoustic tomographic 1" ')       # acoustic material region
#cubit.cmd('block 2 attribute count 2')
#cubit.cmd('block 2 attribute index 1 -1')     # material 1
#cubit.cmd('block 2 attribute index 2 2')      # tomographic model flag


cubit.cmd('export mesh "top.e" dimension 3 overwrite')
cubit.cmd('save as "meshing.cub" overwrite')

#### Export to SESAME format using cubit2specfem3d.py of GEOCUBIT

os.system('mkdir -p MESH')
cubit2specfem3d.export2SESAME('MESH') 

# all files needed by SCOTCH are now in directory MESH
