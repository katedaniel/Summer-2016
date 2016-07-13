import numpy as np
import matplotlib.pyplot as plt

#function that finds root mean square
def rmse(values):
    return np.sqrt((values**2).mean())

dataFilePath = "C:/Users/Noah/Documents/GitHub/Trapped_Orbital_Integrator/table_25.txt"
tableInfo = np.loadtxt(dataFilePath,delimiter=" ",dtype= str)
filepaths = tableInfo[:,0]
data = tableInfo[:,1:16].astype(float)

rms1 = np.sqrt((data[:,[10,11]]**2).mean(axis=1)) #rms 0 to 0.5
rms2 = np.sqrt((data[:,[10,12]]**2).mean(axis=1)) #rms 0 to 1
rms3 = np.sqrt((data[:,[10,13]]**2).mean(axis=1)) #rms 0 to 1.5
rms4 = np.sqrt((data[:,[10,14]]**2).mean(axis=1)) #rms 0 to 2.0

lam_spec = data[:,9]
start_trapped = float(((lam_spec == 0) or (lam_spec == 1) or (lam_spec == 2)).sum())
end_trapped = float(((lam_spec == 0) or (lam_spec == 1)).sum())
trap_frac = end_trapped/start_trapped