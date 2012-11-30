import numpy as np
import os
import sys


num_events_1cable = 3
num_events = num_events_1cable**2

lo = np.linspace(250.0,750.0,num_events_1cable)
la = np.linspace(250.0,750.0,num_events_1cable)
#lo = np.linspace(250.0,2750.0,num_events_1cable)
#la = np.linspace(250.0,2750.0,num_events_1cable)

lon,lat,depth = np.zeros(num_events), np.zeros(num_events), np.zeros(num_events)
for i in range(num_events_1cable):
    for j in range(num_events_1cable):
        lon[num_events_1cable*i+j] = lo[i]
        lat[num_events_1cable*i+j] = la[j]
        # acoustic sources are 0 on the free surface
        # therefore we place them 10 m (0.01 km) deep
        depth[num_events_1cable*i+j] = 0.01 

time_shift = 0.0
halfdur = 0.2

# components of the moment tensor (diag must be equal for acoustic source)
Mrr = 1.0e+23
Mtt = 1.0e+23
Mpp = 1.0e+23
Mrt = 0.0
Mrp = 0.0
Mtp = 0.0

os.system('rm -rf CMTSOLUTION_files')
os.system('mkdir CMTSOLUTION_files')

for i in range(num_events):
    f=open('CMTSOLUTION_files/event'+str(i),'w')    
    f.write('PDE 2011 1 1 1 1 1 '+str(lat[i])+' '+str(lon[i])+' '+str(depth[i])+' 0 0 '+str(i)+'\n')
    f.write('event name:       hom_explosion\n')
    f.write('time shift:       '+str(time_shift)+'\n')
    f.write('half duration:       '+str(halfdur)+'\n')
    f.write('latitute:       '+str(lat[i])+'\n')
    f.write('longitude:       '+str(lon[i])+'\n')
    f.write('depth:       '+str(depth[i])+'\n')
    f.write('Mrr:       '+str(Mrr)+'\n')
    f.write('Mtt:       '+str(Mtt)+'\n')
    f.write('Mpp:       '+str(Mpp)+'\n')
    f.write('Mrt:       '+str(Mrt)+'\n')
    f.write('Mrp:       '+str(Mrp)+'\n')
    f.write('Mtp:       '+str(Mtp)+'\n')
    f.close() # close CMTSOLUTION_files/event#            
        
