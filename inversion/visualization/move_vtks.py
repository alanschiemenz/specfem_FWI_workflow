from glob import glob
import os
vtus = glob('*/*.vtk')

for f in vtus:
    fnew = f[:f.find('/')]+'.vtk'
    os.system('mv '+f+' '+fnew)
