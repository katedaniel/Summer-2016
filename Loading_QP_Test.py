import numpy as np
import os
import Orbit_Code 
reload(Orbit_Code)

#from string of filename, return array of initial conditions
def parseFilename(filename):
    parts = filename.split("=")
    args = np.empty(len(parts)-1)
    for l in range(1, len(parts)):
        num = float(parts[l].split(")")[0])
        args[l-1] = num
    return args
      
#from list of filenames, return matrix of initial conditions    
def parseList(files):
    table = np.empty([len(files),9])
    for l in range(0,len(files)):
        table[l] = parseFilename(files[l])
    return table
             
filepath = "/Users/LBarbano/Desktop/QP_Dump/"
files = os.listdir(filepath)[1:]
table = parseList(files)

a = table[-1,:]

orbit = Orbit_Code.Orbit_Calculator(a[0],a[1],a[2],a[3],a[4],a[5],a[6],a[7],a[8])
f = np.load(filepath + files[-1])

orbit.setqp(f)
orbit.plot(1)
#orbit.Poincare()
 