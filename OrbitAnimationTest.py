import numpy as np
from matplotlib import animation
import Orbit_Code 
reload(Orbit_Code)

def parseFilename(filename):
   
    parts = filename.split("=")
    args = []
    for l in range(1, len(parts)):
        num = float(parts[l].split(")")[0])
        args.append(num)
    return args
    
filepath = "C:/Users/Noah/Documents/GitHub/Trapped_Orbital_Integrator/qp_file/"

#Copy and paste the filename of the desired qp file here
filename = "qp_(m=4)_(th=25.0)_(t=5.0)_(CR=8.0)_(eps=0.4)_(x0=7.8)_(y0=0.0)_(vx0=3.0)_(vy0=230.0).npy"

filename = filepath + filename 
a = parseFilename(filename)
orbit = Orbit_Code.Orbit_Calculator(a[0],a[1],a[2],a[3],a[4],a[5],a[6],a[7],a[8])
f = np.load(filename)
orbit.setqp(f)

fig,ax = orbit.plot(2)

RetardingConstant = 0.001*orbit.getNSteps() #Animation related stuff

star = orbit.getqpR()
guide = np.transpose(orbit.findRg())
#Star and guiding center
s, = ax.plot(star[0,0],star[0,1], 'r*',markersize='12') #Draw initial location of star
g, = ax.plot(guide[0,0],guide[0,1], 'ko',markersize='6') #Draw initial location of star
#Path of star and guiding center
paths, = ax.plot(star[0,0],star[0,1], color="red",ls='dotted')
pathg, = ax.plot(guide[0,0],guide[0,1], color="black",ls='dashed')

# animation function.  This is called sequentially
def animate(i):
    s.set_data(star[int(RetardingConstant*i),0],star[int(RetardingConstant*i),1])  
    paths.set_data(star[0:int(RetardingConstant*i),0],star[0:int(RetardingConstant*i),1])
    
    g.set_data(guide[int(RetardingConstant*i),0],guide[int(RetardingConstant*i),1])  
    pathg.set_data(guide[0:int(RetardingConstant*i),0],guide[0:int(RetardingConstant*i),1])
    
    return s,g,paths,pathg
    
# call the animator.  blit=True means only re-draw the parts that have changed.
anim = animation.FuncAnimation(fig, animate, frames= int(orbit.getNSteps()/RetardingConstant), interval= 1.0)

