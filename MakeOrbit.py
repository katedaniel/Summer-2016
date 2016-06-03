import astropy.units as u
import astropy.constants as const
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from timeit import default_timer
start = default_timer()

#########################
# Defining some constants
#########################
pi = np.pi
G = const.G

###################################
# Defining model dependent constants
###################################

### Disk parameters
# You can toggle these
vc = 220 *u.km /u.s # Circular velocity of the disk
RSun = 8 *u.kpc # Solar galactocentric radius
SigmaSun = 50 *u.Msun /u.pc**2 # Surface density at the solar radius
Rd = 2.5 *u.kpc # Scale length of the disk


# Functions that define disk parameters
def findSurfaceDensity(R): # Calculates the surface density of the disk at any radius
    Sigma0 = SigmaSun *np.exp(RSun/Rd)
    SigmaR = Sigma0 *np.exp(-R/Rd)
    return SigmaR

### Some time dependent things you can toggle
IntTime = 0.5 *u.Gyr # Total time over which to integrate
StepTime = 100000. *u.yr # Time between each calculation
# DON'T TOUCH the next three lines. They make the time series

IntTimeUnitless = (IntTime /u.yr).decompose()
NSteps = np.rint( (IntTime / StepTime).decompose() ) # Intiger number of total steps to take
T = np.linspace(0,IntTimeUnitless,NSteps)

RetardingConstant = 0.0001*NSteps #Animation related stuff
 
### Spiral paremeters 
# You can toggle these  
m = 4
theta = 25. *u.degree
CR = 8. *u.kpc
epsilon = 0.4
OmegaCR = vc/CR # Define orbital frequency



# Functions that define spiral parameters
def findalpha(m,theta): # Calculated parameter alpha
    alpha = m/np.tan(theta)
    return alpha

def findA(x,y,epsilon,alpha): # Calculated amplitude of spiral perturbation
    R = np.sqrt(x**2 + y**2)
    SigmaCR = findSurfaceDensity(R)
    A = 2. *pi *G *SigmaCR*epsilon*R/alpha
    return A

# Use functions to define spiral parameters
alpha = findalpha(m,theta)

#############
# Orbital Integrator
#############

def dvdt(qp,tnow): # Finds the accelleration in this potential at coordinate x-y
    x = qp[0] *u.kpc
    y = qp[1] *u.kpc
    A = findA(x,y,epsilon,alpha)
    time = tnow *u.yr # Note: If you change units, be sure to also change in makeorbit 
    # Find acceleration from logarithmic disk
    dvxD = vc**2 *x/(x**2 + y**2)
    dvyD = vc**2 *y/(x**2 + y**2)
    # Find acceleration from spiral
    dvSFront = -A *np.sin((m *time *vc /CR)*u.rad -m *np.arctan(y/x) -alpha *np.log( np.sqrt( x**2+y**2)/CR)*u.rad)/(x**2 +y**2)
    dvxS = dvSFront *(m*y + (-alpha-(1-np.sqrt(x**2+y**2)/Rd)*(np.tan((m *time *vc /CR)*u.rad -m *np.arctan(y/x) -alpha *np.log( np.sqrt( x**2+y**2)/CR)*u.rad)**-1))*x)
    dvyS = dvSFront *(-m*x + (-alpha-(1-np.sqrt(x**2+y**2)/Rd)*(np.tan((m *time *vc /CR)*u.rad -m *np.arctan(y/x) -alpha *np.log( np.sqrt( x**2+y**2)/CR)*u.rad)**-1))*y)
    # Find total acceleration
    dvxdt = ((dvxD + dvxS)/(u.km /u.s**2))
    dvydt = ((dvyD + dvyS)/(u.km /u.s**2))
    return np.array([dvxdt,dvydt])

def leapstep(qpnow,tnow): # A single leapstep (t+dt), using kick-drift-kick method
    dt = (StepTime/u.year)*3.15576e+07 #convert StepTime to seconds and remove units
    x = qpnow[0]
    y = qpnow[1]
    vx = qpnow[2]
    vy = qpnow[3]
    # x and y are in kpc, vx and vy are in km/s
    a = dvdt(qpnow,tnow) # Find acceleration at this coordinate
    vx = vx -0.5 *dt *a[0] # Advance v_x by half step
    vy = vy -0.5 *dt *a[1] # Advance v_y by half step
    x = x +(dt*vx*3.24077928947e-17) # Advance x by full step, while converting v*dt from km to kpc
    y = y +(dt*vy*3.24077928947e-17) # Advance y by full step, while converting v*dt from km to kpc
    # x and y are in kpc, vx and vy are in km/s
    qpmid = np.array([x,y,vx,vy])
    a  =dvdt(qpmid,tnow) # Find a at new position and complete the velocity step
    vx = vx -0.5 *dt *a[0] # Complete v_x step
    vy = vy -0.5 *dt *a[1] # Complete v_y step
    qpnew = np.array([x,y,vx,vy,tnow]) 
    return qpnew
    
### Finding various theoretical properties

def findEj(qp_part):
    #pulling info
    x = qp_part[:,0]*u.kpc
    y = qp_part[:,1]*u.kpc
    vx = qp_part[:,2]*u.km/u.s
    vy = qp_part[:,3]*u.km/u.s
    t = qp_part[:,4]*u.yr
    R = np.sqrt((x**2)+(y**2))
    phi = np.arctan2(y,x)
    #finding disk potential
    disk_potential = (vc**2)*np.log(R/u.kpc)
    #finding spiral potential
    SigmaR = findSurfaceDensity(R)
    A_r = 2. *pi *G *SigmaR*epsilon*R/alpha
    spiral_potential = A_r.to((u.km/u.s)**2)*np.cos(-alpha*np.log(R/CR)*u.rad + (m*vc*t/CR)*u.rad -m*phi)
    #calculating potential, then energy, then Ej
    potential = disk_potential + spiral_potential
    E_tot = potential + 0.5*(vx**2 + vy**2)
    L_z = R*-np.sqrt(vx**2 + vy**2)*np.sin(phi - np.arctan2(vy,vx))
    E_j = E_tot - OmegaCR*L_z
    return np.array([E_j, L_z, E_tot]) #implicit units of (km/s)**2


def makeorbit(qp0):
    qp = np.zeros(shape=(len(T),5))
    qp[0] = np.append(qp0,[0])
    print NSteps
    for i in range(len(T)):
        qpstep = leapstep(qp0,T[i])
        qp[i] = qpstep
        qp0 = qpstep
    return qp


def toRframe(qp):  # Convert coordinates from NR-frame to R-frame
    # Pull out cartesian non-rotating info
    x = qp[:,0]
    y = qp[:,1]
    vx = qp[:,2]
    vy = qp[:,3]
    t = T *u.yr
    # Calculate polar rotating coordinates for position
    R = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y,x)*u.rad
    phiR= phi - (t*OmegaCR).decompose() *u.rad
    xR = R *np.cos(phiR)
    yR = R *np.sin(phiR)
    # Calculate polar rotating coordinates for velocity
    v_tot = np.sqrt(vx**2 + vy**2)
    theta = np.arctan2(vy,vx)*u.rad #angle between velocity vector and x-axis
    alph = phi - theta # angle between position and velocity vectors
    vr = v_tot*np.cos(alph)
    vphi = -v_tot*np.sin(alph)
    vphiR = vphi - ((OmegaCR*R*u.kpc)/(u.km/u.s))
    vxR = vr *np.cos(phiR) + vphiR*np.sin(phiR)
    vyR = vr *np.sin(phiR) + vphiR*np.cos(phiR)
    # Put it all back into a new array for rotating frame
    qpR = np.transpose(np.array([xR,yR,vxR,vyR,vr,vphi,vphiR]))
    return qpR

def plotArms():
    points = 100
    radius = np.zeros(points)
    t = np.linspace(0,np.pi/2)
    
    for i in xrange(0,m):
    
        radius = (CR/u.kpc)*np.exp((-m*(t) +np.pi)/alpha)
        ax.plot(radius*np.cos(t+2*np.pi*i/m),radius*np.sin(t+2*np.pi*i/m), color="purple",ls='dashed', lw='2.')
  
    return None

####Lindlbad Resonance stuff
#for the first resonances
R_1o = (m+np.sqrt(2))*vc/(m*OmegaCR)
R_1i = (m-np.sqrt(2))*vc/(m*OmegaCR)
#for the ultraharmonic resonances
R_2o = ((2*m)+np.sqrt(2))*vc/((2*m)*OmegaCR)
R_2i = ((2*m)-np.sqrt(2))*vc/((2*m)*OmegaCR)


x0 = 7.8 # Must use implicit units of kpc
y0 = 0. # Must use implicit units of kpc
vx0 = 3. # Must use implicit units of km/s
vy0 = 230. # Must use implicit units of km/s
qp0 = np.array([x0,y0,vx0,vy0])
qp = makeorbit(qp0)
qpR = toRframe(qp)


plt.close('all') #close old plots still up


fig=plt.figure(1) #setting up the basic figure with axes and labels
ax=fig.add_subplot(1,1,1)
plt.xlabel(r'$x$ (kpc)')
plt.ylabel(r'$y$ (kpc)')
plt.axis([-13,13,-13,13])

circ = plt.Circle((0,0), (CR/u.kpc), color='g', fill=False) #plotting CR radius
ax.add_patch(circ)

lind1 = plt.Circle((0,0), (R_1o/u.kpc), color='g', fill=False,ls = 'dashed')
lind2 = plt.Circle((0,0), (R_1i/u.kpc), color='g', fill=False,ls = 'dashed')
lind1uh = plt.Circle((0,0), (R_2o/u.kpc), color='g', fill=False,ls = 'dotted')
lind2uh = plt.Circle((0,0), (R_2i/u.kpc), color='g', fill=False,ls = 'dotted')
ax.add_patch(lind1)
ax.add_patch(lind2)
ax.add_patch(lind1uh)
ax.add_patch(lind2uh)

plt.plot(qpR[:,0][0], qpR[:,1][0], 'g*', markersize='9') #plotting the start of the stellar path
plt.plot(qpR[:,0],qpR[:,1], color="SlateBlue", markevery=500, marker='.', ms=8) 
#plotting the stellar path, markers at (markerevery*StepTime) time, e.g. 10^7

plotArms()

'''
#Animation Stuff
###############################################################################
#Star
s, = ax.plot(qpR[0,0],qpR[0,1], 'r*',markersize='12') #Draw initial location of star
#Path
paths, = ax.plot(qpR[0,0],qpR[0,1], color="red",ls='dotted')

# animation function.  This is called sequentially
def animate(i):
    s.set_data(qpR[int(RetardingConstant*i),0],qpR[int(RetardingConstant*i),1])  
    paths.set_data(qpR[0:int(RetardingConstant*i),0],qpR[0:int(RetardingConstant*i),1])
    return s,paths
    
# call the animator.  blit=True means only re-draw the parts that have changed.
anim = animation.FuncAnimation(fig, animate, frames=int(NSteps/RetardingConstant), interval=0.05, blit = True)
###############################################################################
plotArms()
'''

plt.show()

duration = default_timer() - start
print 'time:'
print duration

filename = "qp_(m=%s)_(t=%s)_(CR=%s)_(eps=%s)_(x0=%s)_(y0=%s)_(vx0=%s)_(vy0=%s)" %(str(m),
str(IntTime/u.Gyr),str(CR/u.kpc),str(epsilon),str(x0),str(y0),str(vx0),str(vy0))

np.save('C:\Summer_2016\qp_file\%s' % filename,qp)