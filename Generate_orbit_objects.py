import numpy as np
from timeit import default_timer
import Orbit_Code as OC
reload(OC)
    
#load the table
theta = 'theta=15'
dataFilePath = "/Users/LBarbano/Github/Summer-2016/table1_"+theta+".txt"
tableInfo = np.loadtxt(dataFilePath,delimiter=" ",dtype= str)
filepaths = tableInfo[:,0]
data = tableInfo[:,1:16].astype(float)

frac= 1 #fraction of qp data to be animated, 1 for all 0.01 for 1%
length = int(frac*len(data))
print("Generating %i Orbit_Code objects for plotting..." % length)
start = default_timer()


orbits = np.array([OC.Orbit_Calculator(data[i,0],data[i,1],data[i,2],data[i,3],data[i,4],
data[i,5],data[i,6],data[i,7],data[i,8]) for i in range(length)]) #Create array of orbit objects
qpdata = np.array([np.loadtxt(filepaths[i],delimiter=" ",dtype= str).astype(float) for i in range(length)]) #extract qpdata from all files
t = qpdata[:,:,0] #switch around t values to so that we can setqp 
qpdata = np.array([np.c_[qpdata[i,:,1:5],t[i,:]] for i in range(length)])

#set qp
for i in range(length):
    orbits[i].setqp(qpdata[i])

qpRdata = np.array([orbits[i].getqpR() for i in range(length)])

duration = default_timer() - start
print "time: %s s" % str(duration) 

#Add time values to qpR
qpR = np.array([np.c_[qpRdata[i,:,0:2],t[i].transpose()] for i in range(length)])

#Save qpR data
np.save('/Users/LBarbano/Github/Summer-2016/qpRdata_(theta=15).npy',qpR)