import numpy as np
import matplotlib.pyplot as plt
from glob import glob
import array
import sys
import pprocess
import os

def read_bin_file(filename, data_types):
    num_objects = len(data_types)
    f = open(filename,'rb')
    file_out=list()
    for n in range(num_objects):
        binlength = array.array('i')
        binlength.read(f,1)
        binvalues = array.array(data_types[n])
        if (data_types[n]=='d'):
            binvalues.read(f,binlength[0]/8) # double precision
        else:
            binvalues.read(f,binlength[0]/4) # 4 for single precision
        binlength.read(f,1)
        file_out.append(binvalues)
    f.close()
    return file_out

def get_plotpoints(p):
    x=read_bin_file(xfiles[p],['f'])[0]
    x=np.array(x,dtype='f')
    y=read_bin_file(yfiles[p],['f'])[0]
    y=np.array(y,dtype='f')
    if not(min(x)-tol<x_coord<max(x)+tol) or not(min(y)-tol<y_coord<max(y)+tol):
        return []
    I=[]
    for gp in range(len(x)):
        if abs(x[gp]-x_coord)<tol and abs(y[gp]-y_coord)<tol:
            I.append(gp)
    z=read_bin_file(zfiles[p],['f'])[0]
    z=np.array(z,dtype='f')
    model=read_bin_file(modelfiles[p],['f'])[0]
    model=np.array(model,dtype='f')
    return model[I],z[I],x[I],y[I]

modeltype='vp' # e.g. vp, alpha_acoustic_kernel, etc

# (x,y) coordinates to center well log on
x_coord = 500.0
y_coord = 500.0

# threshold to consider when searching for (x,y) points near (x_coord,y_coord)
tol=1.0

db_dir='../OUTPUT_FILES/DATABASES_MPI/'
xfiles=glob(db_dir+'*xlocal*')
yfiles=glob(db_dir+'*ylocal*')
zfiles=glob(db_dir+'*zlocal*')
xfiles.sort()
yfiles.sort()
zfiles.sort()
NPROC = len(xfiles)

ll=[]
models_toplot=['true_model']
for mm in [0,3,6,9,10]:
    models_toplot.append('m'+str(mm))

for model_dir in models_toplot:
    modelfiles = glob(model_dir+'/*'+modeltype+'.bin')
    modelfiles.sort()
    results = pprocess.Map(limit=4,reuse=1)
    parallel_function = results.manage(pprocess.MakeReusable(get_plotpoints))
    for procnum in range(NPROC):
        parallel_function(procnum)            
    results.finish()
    
    modelplot,zplot = [],[]
    for procnum in range(NPROC):
        if results[procnum]==[]:
            continue
        else:
            modelplot=np.concatenate((modelplot,results[procnum][0]))
            zplot=np.concatenate((zplot,results[procnum][1]))
    Isort=np.argsort(zplot)
    plt.plot(modelplot[Isort],zplot[Isort],'.-')
    ll.append(model_dir)
    print "Model " + model_dir + " -- finished"
plt.legend(ll)
plt.xlabel('vp (m/s)')
plt.ylabel('z (m)')
print "Saving to file well_log.png"
plt.savefig('well_log.png')
plt.show()
    
    
 
