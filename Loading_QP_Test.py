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
  
filepath = "C:/Trapped_Orbital_Integrator/qp_file/"

#Copy and paste the filename of the desired qp file here

filename = "qp_(m=4)_(th=25.0)_(t=5.0)_(CR=8.0)_(eps=0.4)_(x0=7.8)_(y0=0.0)_(vx0=3.0)_(vy0=230.0).npy"

filename = filepath + filename 
a = parseFilename(filename)
orbit = Orbit_Code.Orbit_Calculator(a[0],a[1],a[2],a[3],a[4],a[5],a[6],a[7],a[8])

f = np.loadtxt(filename)

#orbit.setqp(f)
#orbit.plot(1)
#orbit.Poincare()