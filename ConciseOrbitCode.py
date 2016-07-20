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

###If you want to create a new qp
#make the orbit
orbit = Orbit_Code.Orbit_Calculator(4,20,.3,8,.3,5.,5.,50.,230.) 
orbit.makeOrbit()
#save the orbit
filename = "C:/Users/Noah/Documents/GitHub/Trapped_Orbital_Integrator/qp_file/"
orbit.saveData(filename)

'''
##If you want to load a previous qp
#Copy and paste the filename of the desired qp file here
filename = 'C:/Trapped_Orbital_Integrator/tar_master_20/qp_file_154/qp_(m=4)_(th=20)_(t=2)_(CR=8)_(eps=0.3)_(x0=7.8)_(y0=-3.31767)_(vx0=168.601)_(vy0=-156.033).txt'
#loading and setting the qp
a = parseFilename(filename)
orbit = Orbit_Code.Orbit_Calculator(a[0],a[1],a[2],a[3],a[4],a[5],a[6],a[7],a[8])
f = np.loadtxt(filename)
orbit.setqp(f)
'''