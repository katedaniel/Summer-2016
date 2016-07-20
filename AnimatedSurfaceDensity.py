import numpy as np
from timeit import default_timer
import matplotlib.animation as animation
import matplotlib.pyplot as plt
import Orbit_Code as OC
reload(OC)
    
#load the table
dataFilePath = "/Users/LBarbano/Github/Summer-2016/table1_theta=20.txt"
tableInfo = np.loadtxt(dataFilePath,delimiter=" ",dtype= str)
filepaths = tableInfo[:,0]
data = tableInfo[:,1:16].astype(float)

frac= 0.01 #fraction of qp data to be animated, 1 for all 0.01 for 1%
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

#plots an animation of all qp data
def animateGalaxy():
    plt.close('all')  
    RetardingConstant = 0.001*len(qpRdata[0]) #Animation related stuff
    fig,ax = orbits[0].plot(2)
    s, = ax.plot(qpRdata[:,0,0],qpRdata[:,0,1], 'b*',markersize='5',markeredgecolor = 'black',alpha = 0.5) #Draw initial location of star   
    anim = animation.FuncAnimation(fig, animateScatter, frames= int(len(qpRdata[0])/RetardingConstant), interval= 1.0,fargs = (fig,ax,RetardingConstant,s))
    ax.set_title("Time: %f" % float(t[0,4]/1000000000.0))
    return anim 

#plots an animation of surface density of all qp data
def animateSurfaceDensity(): 
    plt.close('all')  
    RetardingConstant = 0.001*len(qpRdata[0]) #Animation related stuff
    fig,ax = orbits[0].plot(2)
    H, xedges, yedges = np.histogram2d(qpRdata[:,0,0],qpRdata[:,0,1], range=[[-20.,20.0], [-20.,20.0]], bins=(45,45))
    xidx = np.clip(np.digitize(qpRdata[:,0,0], xedges), 0, H.shape[0]-1)
    yidx = np.clip(np.digitize(qpRdata[:,0,1], yedges), 0, H.shape[1]-1)
    col = H[xidx, yidx]
    scat = ax.scatter(qpRdata[:,0,0],qpRdata[:,0,1],c = col, s = 40,cmap=plt.cm.PuBuGn, edgecolor = 'none') #Draw initial location of star 
    ax.set_title("Time: %f" % float(t[0,0]/1000000000.0))
    plt.colorbar(scat)
    anim = animation.FuncAnimation(fig, animateHist, frames= int(len(qpRdata[0])/RetardingConstant), interval= 1.0, 
    fargs=(fig,ax,scat,RetardingConstant,xedges,yedges))
    return anim

#Helper functions for animations, don't call these directly   
def animateScatter(i,fig,ax,RetardingConstant,s):
    s.set_data(qpRdata[:,int(RetardingConstant*i),0],qpRdata[:,int(RetardingConstant*i),1]) 
    ax.set_title("Time(Gy): %f" % float(t[0,int(RetardingConstant*i)]/1000000000.0))
    return s 
         
def animateHist(i,fig,ax,scat,RetardingConstant,xedges,yedges):
    x = qpRdata[:,int(RetardingConstant*i),0]
    y = qpRdata[:,int(RetardingConstant*i),1]
    data = np.hstack((x[:,np.newaxis], y[:, np.newaxis]))
    scat.set_offsets(data)    
    H, xxedges, yyedges = np.histogram2d(qpRdata[:,int(RetardingConstant*i),0],qpRdata[:,int(RetardingConstant*i),1], range=[[-20.,20.0], [-20.,20.0]], bins=(45, 45))
    xidx = np.clip(np.digitize(qpRdata[:,int(RetardingConstant*i),0], xedges), 0, H.shape[0]-1)
    yidx = np.clip(np.digitize(qpRdata[:,int(RetardingConstant*i),1], yedges), 0, H.shape[1]-1)
    col = H[xidx, yidx]
    scat.set_array(col)
    ax.set_title("Time(Gy): %f" % float(t[0,int(RetardingConstant*i)]/1000000000.0))
    return scat


