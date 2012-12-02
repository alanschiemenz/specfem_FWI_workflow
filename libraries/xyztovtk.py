import numpy as np
filename = 'tomography_model.xyz'

ftop = open(filename,'r')
dummy = ftop.readline() # x,y,z min; x,y,z max
dummy = ftop.readline() # dx dy dz
numdim = ftop.readline() # num elements x,y,z
numdim = numdim.split()
numpoints = int(numdim[0])*int(numdim[1])*int(numdim[2])
dummy = ftop.readline() # min/max vp, vs, rho
ftop.close()


f = open('tomography_model.vtk','w')
f.write('# vtk DataFile Version 3.8.1\n')
f.write('material model VTK file\n')
f.write('ASCII\n')
f.write('DATASET STRUCTURED_GRID\n')
f.write('DIMENSIONS '+numdim[0]+' '+numdim[1]+' '+numdim[2]+'\n')
f.write('POINTS '+str(numpoints)+' float\n')

fin = np.loadtxt(filename,skiprows=4,usecols=(0,1,2,3)) # 0,1,2 are x,y,z; 3 is vp
#f.write(fin)

for p in range(len(fin)):
    f.write('%f %f %f \n' % (fin[p,0],fin[p,1],fin[p,2]))

f.write('\n')
f.write('POINT_DATA '+str(numpoints)+'\n')
f.write('SCALARS Vp float 1\n')
f.write('LOOKUP_TABLE default\n')

for p in range(len(fin)):
    f.write('%f \n' % (fin[p,3]))

f.close()
