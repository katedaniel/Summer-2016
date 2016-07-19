import Table_Helper as tableMaker
reload(tableMaker)
import numpy as np
from timeit import default_timer

start = default_timer()
print("Calculating...")

filepath = "C:/Trapped_Orbital_Integrator/tar_master_30/" #filepath to tar files
#tableMaker.unTar(filepath) #untar the files
table, table2 = tableMaker.genTable(filepath) #make table of analysis stuff

duration = default_timer() - start
print "time: %s s" % str(duration) 
'''
#This table is for individual qp's with some of their data
dataFilePath = "C:/Users/Noah/Documents/GitHub/Trapped_Orbital_Integrator/table_30.txt"
np.savetxt(dataFilePath, table.astype(str), delimiter = " ",fmt='%s')
'''
#This table is for accumulated data on trapped fraction and angmom rms
dataFilePath2 = "C:/Users/Noah/Documents/GitHub/Trapped_Orbital_Integrator/table2_30.txt"
np.savetxt(dataFilePath2, table2.astype(str), delimiter = " ",fmt='%s')