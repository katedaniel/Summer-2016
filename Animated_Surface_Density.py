import numpy as np
from timeit import default_timer
import matplotlib.animation as animation
import matplotlib.pyplot as plt
import Orbit_Code as OC
reload(OC)
    
print("Load stuff for plotting...")
start = default_timer()

qpRdata = np.load('/Users/LBarbano/Desktop/QP_Data/qpRdata_(theta=20).npy')
#Create orbit object to plot the spiral parameters
#Hardcode spiral parameters (m=4,theta=15,20,25,30 etc, t = ?, eps = 0.3, rest don't matter) 
orbit = OC.Orbit_Calculator(4,30,2,8,0.3,0,0,0,0) 
duration = default_timer() - start
print "time: %s s" % str(duration) 


#plots an animation of all qp data
def animateGalaxy():
    plt.close('all')  
    RetardingConstant = 0.001*len(qpRdata[0]) #Animation related stuff
    fig,ax = orbit.plot(2)
    s, = ax.plot(qpRdata[:,0,0],qpRdata[:,0,1], 'b*',markersize='3.5',markeredgecolor = 'black',alpha = 1) #Draw initial location of star   
    anim = animation.FuncAnimation(fig, animateScatter, frames= int(len(qpRdata[0])/RetardingConstant), interval= 1.0,fargs = (fig,ax,RetardingConstant,s))
    ax.set_title("Time: %f" % float(qpRdata[0,0,2]/1000000000.0))
    return anim 

#plots an animation of surface density of all qp data
def animateSurfaceDensity(): 
    plt.close('all')  
    RetardingConstant = 0.001*len(qpRdata[0]) #Animation related stuff
    fig,ax = orbit.plot(2)
    H, xedges, yedges = np.histogram2d(qpRdata[:,0,0],qpRdata[:,0,1], range=[[-15.,15.0], [-15.,15.0]], bins=(40,40))
    xidx = np.clip(np.digitize(qpRdata[:,0,0], xedges+.75), 0, H.shape[0]-1)
    yidx = np.clip(np.digitize(qpRdata[:,0,1], yedges+.75), 0, H.shape[1]-1)
    col = H[xidx, yidx]
    scat = ax.scatter(qpRdata[:,0,0],qpRdata[:,0,1],c = col, s = 8,cmap='jet', edgecolor = 'none') #Draw initial location of star 
    ax.set_title("Time: %f" % float(qpRdata[0,0,2]/1000000000.0))
    plt.colorbar(scat)
    anim = animation.FuncAnimation(fig, animateHist, frames= int(len(qpRdata[0])/RetardingConstant), interval= 1.0, 
    fargs=(fig,ax,scat,RetardingConstant,xedges,yedges))
    return anim

#Helper functions for animations, don't call these directly   
def animateScatter(i,fig,ax,RetardingConstant,s):
    s.set_data(qpRdata[:,int(RetardingConstant*i),0],qpRdata[:,int(RetardingConstant*i),1]) 
    ax.set_title("Time(Gy): %f" % float(qpRdata[0,int(RetardingConstant*i),2]/1000000000.0))
    return s 
         
def animateHist(i,fig,ax,scat,RetardingConstant,xedges,yedges):
    x = qpRdata[:,int(RetardingConstant*i),0]
    y = qpRdata[:,int(RetardingConstant*i),1]
    data = np.hstack((x[:,np.newaxis], y[:, np.newaxis]))
    scat.set_offsets(data)    
    H, xxedges, yyedges = np.histogram2d(qpRdata[:,int(RetardingConstant*i),0],qpRdata[:,int(RetardingConstant*i),1], range=[[-15.,15.0], [-15.,15.0]], bins=(40, 40))
    xidx = np.clip(np.digitize(qpRdata[:,int(RetardingConstant*i),0], xedges+.75), 0, H.shape[0]-1)
    yidx = np.clip(np.digitize(qpRdata[:,int(RetardingConstant*i),1], yedges+.75), 0, H.shape[1]-1)
    col = H[xidx, yidx]
    scat.set_array(col)
    ax.set_title("Time(Gy): %f" % float(qpRdata[0,int(RetardingConstant*i),2]/1000000000.0))
    return scat


def plotSurfaceDensity(i):
  
    fig,ax = orbit.plot(2)
    H, xedges, yedges = np.histogram2d(qpRdata[:,i,0],qpRdata[:,i,1], range=[[-15.,15.0], [-15.,15.0]], bins=(40,40))
    xidx = np.clip(np.digitize(qpRdata[:,i,0], xedges+0.75), 0, H.shape[0]-1)
    yidx = np.clip(np.digitize(qpRdata[:,i,1], yedges+0.75), 0, H.shape[1]-1)
    col = H[xidx, yidx]
    scat = ax.scatter(qpRdata[:,i,0],qpRdata[:,i,1],c = col, s = 8,cmap='jet', edgecolor = 'none') #Draw initial location of star 
    ax.set_title("Time(Gy): %f" % float(qpRdata[0,i,2]/1000000000.0))
    plt.colorbar(scat)
    plt.show()