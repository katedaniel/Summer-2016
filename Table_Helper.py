import numpy as np
import os
import Orbit_Code 
reload(Orbit_Code)
import tarfile

#from string of filename, return array of initial conditions
def parseFilename(filename):
    parts = filename.split("=")
    args = np.empty(len(parts)-1)
    for l in range(1, len(parts)):
        num = float(parts[l].split(")")[0])
        args[l-1] = num
    return args
      
#from list of filenames, return matrix of initial conditions    
def parseList(files):
    table = np.empty([len(files),9])
    for l in range(0,len(files)):
        table[l] = parseFilename(files[l])
    return table
 
def genTable(filepath):
    table = []
    qp_len = 2001.       #the length of qp, can't figure out how to code that in
    table2 = np.zeros(qp_len)  
    table3 = np.zeros(qp_len)
    table4 = np.zeros(qp_len)  
    for dirpath, dirnames, files in os.walk(filepath):
        for f in files:
            fullpath = os.path.join(dirpath, f) #get full path of subject file
            if f != ".DS_Store" and not f.endswith("tar.gz"):
                a = parseFilename(f) #Turn filename into numpy array of initial conditions
                orbit = Orbit_Code.Orbit_Calculator(a[0],a[1],a[2],a[3],a[4],a[5],a[6],a[7],a[8])
                fullpath = os.path.join(dirpath, f) #get full path of subject file
                data = np.loadtxt(fullpath) #need to change around order of data columns for real thing
                data = data.astype(float)   #change to float
                t = data[:,0]   #next two lines switch order of t,x,y,vx,vy to x,y,vx,vy,t
                data = np.c_[data[:,1:5] ,t] 
                data[:,2] = data[:,2]*9.777922216731282e+8
                data[:,3] = data[:,3]*9.777922216731282e+8
                orbit.setqp(data)
                lam = orbit.findLam()[0]
                if np.absolute(lam[0]) < 1.:
                       table2 += (np.absolute(lam) < 1.)
                       angmom = orbit.findLam()[3]
                       angmom_del = angmom - angmom[0]
                       angmom_del = angmom_del**2
                       table3 += angmom_del
                       if ((np.absolute(lam) < 1).sum() == len(lam)):
                           table4 +=  angmom_del                    
                lamsp = orbit.Lam_special()
                Lz = orbit.findLz()
                table.append([dirpath+'/'+f,a[0],a[1],a[2],a[3],a[4],a[5],a[6],a[7],a[8],lamsp,Lz[0],Lz[1],Lz[2],Lz[3],Lz[4]]) 
    table3 = np.sqrt(table3/qp_len)
    table4 = np.sqrt(table4/qp_len)
    table2 = table2/(table2[0])
    table_final = np.vstack((table3,table4))
    table_final = np.vstack((table2,table_final))
    table_final = np.vstack((t,table_final))
    return np.array(table), table_final.transpose()

def unTar(filepath):
    for dirpath, dirnames, files in os.walk(filepath):
        for f in files:
            fullpath = os.path.join(dirpath, f) #get full path of subject file
            if f.endswith("tar.gz"):
                tar = tarfile.open(fullpath, "r:gz")
                tar.extractall(path = filepath)
                tar.close()
    return