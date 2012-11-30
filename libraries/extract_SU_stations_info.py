# Author: Alan Schiemenz (schiemenz@geophysik.uni-muenchen.de) 
# input: SU_stations_info.bin
# output: SU_stations_coords_proc.txt
#   This script generates an N x 4 array, where N is the total number of stations.  The 4 columns give x,y, and z-coordinates of the stations, and the processor index i of the SU output, i.e. "<i>_dz_SU"

import array

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

SU_stations_info = read_bin_file('SU_stations_info.bin',['i','d','d','d'])

# SU_stations_info[2] contains all x,y,z coords of stations
nstations = len(SU_stations_info[2])/3 
sx=SU_stations_info[2][:nstations]
sy=SU_stations_info[2][nstations:2*nstations]
sz=SU_stations_info[2][2*nstations:3*nstations]
station_CPU_index = SU_stations_info[0][:nstations]

station_out = open('SU_stations_coords_proc.txt','w')
for n in range(nstations):
    station_out.write('%f %f %f %i \n' % (sx[n],sy[n],sz[n],station_CPU_index[n]))
