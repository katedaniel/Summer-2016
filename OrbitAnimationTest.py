import numpy as np
from matplotlib import animation
from Orbit_Calculator import Orbit_Calculator

orbit = Orbit_Calculator(4,0.01,8,.4,7.8,0.,3.,230.) 
orbit.makeOrbit()
fig,ax = orbit.plot(2)

RetardingConstant = 0.005*orbit.getNSteps() #Animation related stuff


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
anim = animation.FuncAnimation(fig, animate, frames= int(orbit.getNSteps()/RetardingConstant), interval= 0.01)
