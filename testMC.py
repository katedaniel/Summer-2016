import numpy as np
import matplotlib.pyplot as plt

plt.close('all')
dim = 16 #dimensions of scatter plot and interval of radii to consider 

#function that finds root mean square
def rmse(values):
    return np.sqrt((values**2).mean())

def sort (radius,velocity):    
    rmsOut = np.array([0,0])
    count = 1 #increments the while loop
    index = 0 #index of lower bound 
    tempRad = radius[0]    
    #sweep through sorted indices to determine how to group velocities to calc rms
    while(count < len(radius)):
        if radius[count] > tempRad:
            rms = rmse(velocity[index:count]) #calc rms          
            rmsOut = np.vstack((rmsOut,[tempRad,rms])) 
            tempRad =  tempRad + 0.5 
            index = count
        count = count +1
    return rmsOut[2:,:]

#Specify filepath to the initial conditions data
filepath = "/Users/LBarbano/Github/Summer-2016/table1_theta=20.txt"
tableInfo = np.loadtxt(filepath,delimiter=" ",dtype= str)
filepaths = tableInfo[:,0]
data = tableInfo[:,1:16].astype(float)

x = data[:,5]
y = data[:,6]
vx = data[:,7]
vy = data[:,8]

'''
Calculate radii to plot histogram of initial positions
By inspection of the histogram, we will be able to determine the exponential 
decay of the disk's surface desnity
'''
#find radius of each initial condition
radius = np.sqrt(x**2 + y**2)
#Plot the histogram!
fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
bins = np.arange(0,100,0.5)
ax1.hist(radius, bins,range=[3, dim])
ax1.set_xlim([3, dim])
ax1.set_xlabel(r' radius (kpc)')
ax1.set_ylabel(r' #of stars in annulus with 1kpc width ')

'''
Calculate velocity dispersion
'''
phi = np.arctan2(y,x) #angle of position vector
v_tot = np.sqrt(vx**2 + vy**2) #total velocity
theta = np.arctan2(vy,vx) #angle between velocity vector and x-axis
alph = phi - theta # angle between position and velocity vectors
vr = v_tot*np.cos(alph)

sortIndex = np.argsort(radius) #indices of sorted radii
radius = np.sort(radius) #actually sort radii
vr = vr[sortIndex] #rearrange vr to correspond to sorted radii
rms = sort(radius,vr) #find rms of annuli with width 0.1kpc

fig2 = plt.figure() #setting up the basic figure with axes and labels
ax2 = fig2.add_subplot(1,1,1)
ax2.set_xlabel(r'radius (kpc)')
ax2.set_ylabel(r'RMS(Vr)')
ax2.scatter(rms[:,0],rms[:,1])
ax2.set_xlim([4, dim])
#ax2.set_ylim([0, 60])


#Plot the distribution!
fig3 = plt.figure()      #setting up the basic figure with axes and labels
ax3 = fig3.add_subplot(1,1,1)
ax3.scatter(x,y, alpha = 0.5)
ax3.set_xlabel(r'$x$ (kpc)')
ax3.set_ylabel(r'$y$ (kpc)')
ax3.axis([-dim,dim,-dim,dim])
ax3.set_aspect('equal', 'datalim')


fig4 = plt.figure()      #setting up the basic figure with axes and labels
ax4 = fig4.add_subplot(1,1,1)
ax4.set_xlabel(r' radius (kpc)')
ax4.set_ylabel(r'vr (km/s)')
ax4.scatter(radius,vr, alpha = 0.3)
ax4.set_xlim([0, dim])


fig5 = plt.figure()      #setting up the basic figure with axes and labels
ax5 = fig5.add_subplot(1,1,1)
ax5.scatter(phi,theta)

bins = 100
fig6 = plt.figure()      #setting up the basic figure with axes and labels
ax6 = fig6.add_subplot(1,1,1)
ax6.hist(alph,bins)


fig7 = plt.figure()
ax7 = fig7.add_subplot(111)
H, xedges, yedges = np.histogram2d(radius, vr, range=[[0,16.0], [-120.,120.0]], bins=(50, 50))
myextent  = [xedges[0],xedges[-1],yedges[0],yedges[-1]]
im = ax7.imshow(H.T,origin='low',extent=myextent,interpolation='nearest',aspect ='auto')
cbar = plt.colorbar(im)
ax7.set_xlim(0, 15)
plt.show()

fig8 = plt.figure()
ax8 = fig8.add_subplot(111)
H, xedges, yedges = np.histogram2d(x, y, range=[[-16.,16.0], [-16.,16.0]], bins=(80, 80))
myextent  = [xedges[0],xedges[-1],yedges[0],yedges[-1]]
im = ax8.imshow(H.T,origin='low',extent=myextent,interpolation='nearest',aspect='equal')
cbar = plt.colorbar(im)

plt.show()


