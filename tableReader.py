import numpy as np
import matplotlib.pyplot as plt

#function that finds root mean square
def rmse(values):
    return np.sqrt((values**2).mean())

dataFilePath = "/Users/LBarbano/Github/Summer-2016/table_theta=30.txt"
tableInfo = np.loadtxt(dataFilePath,delimiter=" ",dtype= str)
filepaths = tableInfo[:,0]
data = tableInfo[:,1:16].astype(float)

rms1 = np.sqrt((data[:,[10,11]]**2).mean(axis=1)) #rms 0 to 0.5
rms2 = np.sqrt((data[:,[10,12]]**2).mean(axis=1)) #rms 0 to 1.0
rms3 = np.sqrt((data[:,[10,13]]**2).mean(axis=1)) #rms 0 to 1.5
rms4 = np.sqrt((data[:,[10,14]]**2).mean(axis=1)) #rms 0 to 2.0

fracStayTrapped = ((data[:,9]==0).sum()+(data[:,9]==1).sum())/((data[:,9]==0).sum()+(data[:,9]==1).sum()+(data[:,9]==2).sum())

zero = (data[:,9]==0).sum()
one = (data[:,9]==1).sum()
two = (data[:,9]==2).sum()
three = (data[:,9]==3).sum()
four = (data[:,9]==4).sum()
five = (data[:,9]==5).sum()

plt.close('all')
fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
bins = np.arange(0,6,1)
ax1.hist(data[:,9], bins,range=[-1, 6])
plt.show()


