import Table_Helper as tableMaker
reload(tableMaker)
import numpy as np
from timeit import default_timer

start = default_timer()
print("Calculating...")

filepath = "/Users/LBarbano/Desktop/QP_Data/Trapped_Orbital_Integrator-master_(theta=30)" #filepath to tar files
#tableMaker.unTar(filepath) #untar the files
table1, table2 = tableMaker.genTable(filepath) #make table of analysis stuff

duration = default_timer() - start
print "Calculation time: %s s" % str(duration) 

dataFilePath1 = "/Users/LBarbano/Github/Summer-2016/table1_theta=30.txt"
np.savetxt(dataFilePath1, table1.astype(str), delimiter = " ",fmt='%s')  

dataFilePath2 = "/Users/LBarbano/Github/Summer-2016/table2_theta=30.txt"
np.savetxt(dataFilePath2, table2.astype(str), delimiter = " ",fmt='%s') 
