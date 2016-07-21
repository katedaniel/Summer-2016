#SurfaceDensityPlot.py plots the surface density of the galaxy given x0 and y0

import numpy as np
import matplotlib.pyplot as plt
from timeit import default_timer

dataFilePath = "C:/Users/Noah/Documents/GitHub/Trapped_Orbital_Integrator/table1_15.txt"
tableInfo = np.loadtxt(dataFilePath,delimiter=" ",dtype= str)
filepaths = tableInfo[:,0]
data = tableInfo[:,1:16].astype(float)
x = data[:,5]
y = data[:,6]

start = default_timer()
plt.close('all')

fig = plt.figure()
ax = fig.add_subplot(111)
H, xedges, yedges = np.histogram2d(x, y, range=[[-16.,16.0], [-16.,16.0]], bins=(80, 80))
myextent  = [xedges[0],xedges[-1],yedges[0],yedges[-1]]
im = ax.imshow(H.T,origin='low',extent=myextent,interpolation='nearest',aspect='equal')
cbar = plt.colorbar(im)
plt.show()

duration = default_timer() - start
print "time: %s s" % str(duration) 