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
    for dirpath, dirnames, files in os.walk(filepath):
        for f in files:
            fullpath = os.path.join(dirpath, f) #get full path of subject file
            if f != ".DS_Store" and not f.endswith("tar.gz"):
                a = parseFilename(f) #Turn filename into numpy array of initial conditions
                orbit = Orbit_Code.Orbit_Calculator(a[0],a[1],a[2],a[3],a[4],a[5],a[6],a[7],a[8])
                fullpath = os.path.join(dirpath, f) #get full path of subject file
                data = np.loadtxt(fullpath)
                data = data.astype(float)   #change to float
                t = data[:,0]   #next two lines switch order of t,x,y,vx,vy to x,y,vx,vy,t
                data = np.c_[data[:,1:5] ,t] 
                data[:,2] = data[:,2]*9.777922216731282e+8
                data[:,3] = data[:,3]*9.777922216731282e+8
                orbit.setqp(data)  
                lamsp = orbit.Lam_special()
                Lz = orbit.findLz()
                #table elements: filepath,m,theta,IntTime,CR,epsilon,x0,y0,vx0,vy0,lamspecial,Lz(0),Lz(0.5),Lz(1.0),Lz(1.5), Lz(2.0)
                table.append([dirpath+"/"+f,a[0],a[1],a[2],a[3],a[4],a[5],a[6],a[7],a[8],lamsp,Lz[0],Lz[1],Lz[2],Lz[3],Lz[4]]) 
    return np.array(table)

def unTar(filepath):
    for dirpath, dirnames, files in os.walk(filepath):
        for f in files:
            fullpath = os.path.join(dirpath, f) #get full path of subject file
            if f.endswith("tar.gz"):
                tar = tarfile.open(fullpath, "r:gz")
                tar.extractall(path = filepath)
                tar.close()
    return
