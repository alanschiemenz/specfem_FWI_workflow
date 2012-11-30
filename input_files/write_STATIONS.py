import numpy as np
import os

print "Writing STATIONS file: \n"
station_file = open('STATIONS',"w")

# equidistant grid of 100 receivers
xR = np.linspace(150.0,850.0,10)
yR = np.linspace(150.0,850.0,10)
#xR = np.linspace(150.0,2850.0,10)
#yR = np.linspace(150.0,2850.0,10)

# Assuming an OBC experiment in 200-m deep water
depth=200.0 

station_num = 0
for yy in range(len(yR)):
    for xx in range(len(xR)):
        xrec, yrec = xR[xx], yR[yy]
        cablenum = str(yy)
        rec_on_cable = str(xx)
        # important : STATIONS file puts y-coordinate before x-coordinate !
        station_file.write("%s %s %f %f %f %f\n" % ('A'+cablenum, rec_on_cable, yrec, xrec, 0.0, depth)) # stations in acoustic region
        station_num = station_num + 1

station_file.close()
print "Output " + str(station_num) + ' receivers to STATIONS file\n'

