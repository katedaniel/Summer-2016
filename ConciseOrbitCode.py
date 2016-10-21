import numpy as np
import Orbit_Code as OC
reload(OC)

def parseFilename(filename):
   
    parts = filename.split("=")
    args = []
    for l in range(1, len(parts)):
        num = float(parts[l].split(")")[0])
        args.append(num)
    return args

###If you want to create a new qp
#make the orbit
orbit = OC.Orbit_Calculator(4,30,.1,8,.3,7.1,0.,90.,220.) 
orbit.makeOrbit()
#save the orbit
filename = "C:/Users/Noah/Documents/GitHub/Trapped_Orbital_Integrator/qp_file/"
#orbit.saveData(filename)

'''
##If you want to load a previous qp
#Copy and paste the filename of the desired qp file here
filename = 'C:/Trapped_Orbital_Integrator/tar_master_20/qp_file_154/qp_(m=4)_(th=20)_(t=2)_(CR=8)_(eps=0.3)_(x0=7.8)_(y0=-3.31767)_(vx0=168.601)_(vy0=-156.033).txt'
#loading and setting the qp
a = parseFilename(filename)
orbit = OC.Orbit_Calculator(a[0],a[1],a[2],a[3],a[4],a[5],a[6],a[7],a[8])
#for home made qp use these lines:
f = np.loadtxt(filename)
orbit.setqp(f)
#for osg made qp use these lines:
#data = np.loadtxt(filename,delimiter=" ",dtype= str).astype(float)
#t = data[:,0]   #next two lines switch order of t,x,y,vx,vy to x,y,vx,vy,t
#data = np.c_[data[:,1:5] ,t] 
#orbit.setqp(data) 
'''