import numpy as np
import Orbit_Code 
reload(Orbit_Code)

def parseFilename(filename):
   
    parts = filename.split("=")
    args = []
    for l in range(1, len(parts)):
        num = float(parts[l].split(")")[0])
        args.append(num)
    return args
'''   
filepath = "C:/Trapped_Orbital_Integrator/qp_file/"

#Copy and paste the filename of the desired qp file here
filename = "qp_(m=4)_(th=25.0)_(t=0.5)_(CR=8.0)_(eps=0.4)_(x0=6.5)_(y0=1.8)_(vx0=-10.0)_(vy0=223.0).npy"

filename = filepath + filename 
a = parseFilename(filename)
orbit = Orbit_Code.Orbit_Calculator(a[0],a[1],a[2],a[3],a[4],a[5],a[6],a[7],a[8])
f = np.load(filename)
'''

orbit = Orbit_Code.Orbit_Calculator(4,25.,2.,8.,0.3,7.8,0.,3.,223.)
f = open('C:/Trapped_Orbital_Integrator/qp_file/qp_(m=4)_(th=25)_(t=2)_(CR=8)_(eps=0.3)_(x0=7.8)_(y0=0)_(vx0=3)_(vy0=223).dat','r')