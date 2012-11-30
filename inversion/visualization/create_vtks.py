from glob import glob
from numpy import array
import os
import sys
import pprocess

def model_loop(model_dir):
    for filetype in filenames: # e.g. vp, alpha_acoustic_kernel
        files = glob(model_dir+'/proc*_'+filetype+'.bin')
        procnum = len(files)
        if procnum>0:
            #fname = dirname+kernel
            # convert .bin files to .mesh format
            # a '0' in the last line of the following command gives coarse output wherease a '1' gives fine output (i.e. every GLL point is shown)
            root_dir=os.getcwd()
            #os.chdir(dirname)                
            cmd = root_dir+'/../bin/xcombine_vol_data 0 '+ str(procnum-1) + ' ' + filetype +' ' + model_dir + ' ' + model_dir + ' 1'
            print cmd
            os.system(cmd)
            #os.chdir(root_dir)

filenames = ['vp']

model_dirs=['true_model/']
for i in range(20):
    model_dirs.append('m'+str(i)+'/')
results = pprocess.Map(limit=len(model_dirs), reuse=1)
parallel_function = results.manage(pprocess.MakeReusable(model_loop))
for model in model_dirs:
    parallel_function(model)
results.finish()

print "Moving vtk files"
vtus=[]
for d in model_dirs: vtus += glob(d+'*.vtk')
for f in vtus:
    fnew = f[:f.find('/')]+'.vtk'
    os.system('mv '+f+' '+fnew)
