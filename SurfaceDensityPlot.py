import numpy as np
import matplotlib.pyplot as plt


dataFilePath = "/Users/LBarbano/Github/Summer-2016/table_theta=30.txt"
tableInfo = np.loadtxt(dataFilePath,delimiter=" ",dtype= str)
filepaths = tableInfo[:,0]
data = tableInfo[:,1:16].astype(float)
x = data[:,5]
y = data[:,6]

plt.close('all')

fig = plt.figure()
ax = fig.add_subplot(111)
H, xedges, yedges = np.histogram2d(x, y, range=[[-20.,20.0], [-20.,20.0]], bins=(50, 50))
myextent  =[xedges[0],xedges[-1],yedges[0],yedges[-1]]
plt.imshow(H.T,origin='low',extent=myextent,interpolation='nearest',aspect='equal')
plt.colorbar()
plt.show()



