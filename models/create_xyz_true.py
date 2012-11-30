from numpy import zeros
import numpy as np

# create file
d = open("tomography_model_true.xyz", "w")

# origin points
ORIG_X=0.
ORIG_Y=0. 
ORIG_Z=0.
# end points
END_X=1000.
END_Y=1000.
END_Z=-1000.  # depth in negative z-direction

# mid point of domain, used to define perturbation
xmid=0.5*(END_X-ORIG_X)
ymid=0.5*(END_Y-ORIG_Y)
zmid=0.5*(END_Z-ORIG_Z)

d.write(str(ORIG_X) + " " + str(ORIG_Y) + " " + str(ORIG_Z) +
        " " + str(END_X) + " " + str(END_Y) + " " + str(END_Z) + "\n")

dx = 25.
dy = 25.
dz = -25.
d.write(str(dx) + " " + str(dy) + " " + str(dz) + "\n") 
# number of interfaces
NX=int((END_X-ORIG_X)/dx)+1
NY=int((END_Y-ORIG_Y)/dy)+1
NZ=int((END_Z-ORIG_Y)/dz)+1

d.write(str(NX) + " " + str(NY) + " " + str(NZ) + "\n") 

VP_MIN = 1000.0
VP_MAX = 4000.0

# constant-density, acoustic
VS_MIN = 0.0
VS_MAX = 0.0
RHO_MIN= 2000.0
RHO_MAX= 2000.0

d.write(str(VP_MIN) + " " + str(VP_MAX) + " " + str(VS_MIN) +
        " " + str(VS_MAX) + " " + str(RHO_MIN) + " " + str(RHO_MAX) + "\n") 

for i in range(NZ):
    for j in range(NY):
        for k in range(NX):
            x = ORIG_X+k*dx
            y = ORIG_Y+j*dy
            z = ORIG_Z+i*dz
            if ((-200 > z > -400)):
                vp = 1480.0 -0.2*z - 250.0 * np.exp(-2.0e-5*((x-xmid)**2+(y-ymid)**2))
            else:
                vp = 1480.0-0.2*z
            # for vs, rho approximations we could Muyzert (2007), Geophysics, eqs 3 and 4
            #vs = 0.0 #0.8621*vp - 1172.4
            #rho = 310*(vp**0.25)

            # constant-density, acoustic
            vs = 0.0
            rho = 2000.0

            d.write(str(x) + " " + str(y) + " " + str(z) + " " +
                    str(vp) + " " + str(vs) + " " + str(rho) + "\n")

d.close()
