import astropy.units as u
import astropy.constants as const
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
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
def findSurfaceDensity(R): # Calculats the surface density of the disk at any radius
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
 
### Spiral paremeters 
# You can toggle these  
m = 4
theta = 25 *u.degree
CR = 8 *u.kpc
epsilon = 0.4

# Functions that define spiral parameters
def findalpha(m,theta): # Calculated parameter alpha
    alpha = m/np.tan(theta)
    return alpha

def findA(CR,epsilon,alpha): # Calculated amplitude of spiral perturbation
    SigmaCR = findSurfaceDensity(CR)
    A = 2. *pi *G *SigmaCR*epsilon*CR/alpha
    return A

# Use functions to define spiral parameters
alpha = findalpha(m,theta)
A = findA(CR,epsilon,alpha)

#############
# Orbital Integrator
#############

def dvdt(qp,tnow): # Finds the accelleration in this potential at coordinate x-y
    x = qp[0] *u.kpc
    y = qp[1] *u.kpc
    time = tnow *u.yr # Note: If you change units, be sure to also change in makeorbit 
    # Find acceleration from logarithmic disk
    dvxD = vc**2 *x/(x**2 + y**2)
    dvyD = vc**2 *y/(x**2 + y**2)
    # Find acceleration from spiral
    dvSFront = -A *np.sin(m *time *vc /CR*u.rad -m *np.arctan(y/x) -alpha *np.log( np.sqrt( x**2+y**2)/CR)*u.rad)/(x**2 +y**2)
    dvxS = dvSFront *(m*y - alpha*x)
    dvyS = -dvSFront *(m*x + alpha*y)
    # Find total acceleration
    dvxdt = ( (dvxD + dvxS)/(u.km /u.s**2) ).decompose()
    dvydt = ( (dvyD + dvyS)/(u.km /u.s**2) ).decompose()
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
    qpnew = np.array([x,y,vx,vy]) 
    return qpnew
    
def makeorbit(qp0):
    qp = np.array([qp0]) #unnecessary, its already an array
    print NSteps
    for x in T:
        qpstep = leapstep(qp0,x)
        qp = np.concatenate((qp,np.array([qpstep])),axis=0)
        qp0 = qpstep
    qp = np.delete(qp,np.shape(qp)[0]-1,0) # Ensure the shape matches for later fcns
    return qp

def findphiR(x,y,t):
    OmegaCR = vc/CR # Define orbital frequency
    rotationangle = (t*OmegaCR).decompose() *u.rad
    if (x < 0  and y >= 0): phiR = np.arctan(y/x) + pi *u.rad - rotationangle
    elif (x < 0  and y < 0):  phiR = np.arctan(y/x) - pi *u.rad - rotationangle
    elif (x == 0 and y > 0):  phiR = pi/2. *u.rad - rotationangle
    elif (x == 0 and y < 0):  phiR = -pi *u.rad - rotationangle
    elif (x == 0 and y == 0): phiR = - rotationangle
    else: phiR = np.arctan(y/x) - rotationangle
    return phiR


def toRframe(qp):  # Convert coordinates from N-frame to R-frame
    OmegaCR = vc/CR # Define orbital frequency
    # Transform into cylindrical coordinates
    x = qp[:,0] *u.kpc
    y = qp[:,1] *u.kpc
    vx = qp[:,2] *u.km/u.s
    vy = qp[:,3] *u.km/u.s
    t = T *u.yr
    R = np.sqrt(x**2 + y**2)
    phiR = np.arctan(y/x) - (t*OmegaCR).decompose() *u.rad #this is a slightly unecessary way to get the array shape, just use qp?
    for i in xrange(0,np.shape(phiR)[0]):
        phiR[i] = findphiR(x[i],y[i],t[i])
    xR = R *np.cos(phiR)
    yR = R *np.sin(phiR)
    vphiR = np.abs(vy *np.cos(np.arctan(y/x))) + np.abs(vx *np.sin(np.arctan(y/x)))
    vxR = vphiR *np.cos(phiR)
    vyR = vphiR *np.sin(phiR)

    xR = (xR/u.kpc).decompose()
    yR = (yR/u.kpc).decompose()
    vxR = (vxR /u.km*u.s).decompose()
    vyR = (vyR /u.km*u.s).decompose()
    R = (R/u.kpc).decompose()
    vphiR = (vphiR/u.km*u.s).decompose()
    qpR = np.transpose(np.array([xR,yR,vxR,vyR,R,vphiR]))
    return qpR



x0 = 7.6 # Must use implicit units of kpc
y0 = 0. # Must use implicit units of kpc
vx0 = -2.5 # Must use implicit units of km/s
vy0 = 223. # Must use implicit units of km/s
qp0 = np.array([x0,y0,vx0,vy0])
qp = makeorbit(qp0)
qpR = toRframe(qp)

'''
plt.xlabel(r'$x$ (kpc)')
plt.ylabel(r'$y$ (kpc)')
plt.axis([-10,10,-10,10])
plt.plot(qp[:,0],qp[:,1], color="SlateBlue")
plt.show()
'''

plt.close('all') #close old plots still up

fig=plt.figure(1) #setting up the basic figure with axes and labels
ax=fig.add_subplot(1,1,1)
plt.xlabel(r'$x$ (kpc)')
plt.ylabel(r'$y$ (kpc)')
plt.axis([-10,10,-10,10])

plt.plot(qpR[:,0],qpR[:,1], color="SlateBlue") #plotting the stellar path
plt.plot(qpR[:,0][0], qpR[:,1][0], 'g*', markersize='12') #plotting the start of the stellar path
circ = plt.Circle((0,0), (CR/u.kpc), color='g', fill=False) #plotting CR radius
linblad1 = plt.Circle((0,0), (CR/u.kpc)*0.8, color='r', fill=False, ls='dashed')
linblad2 = plt.Circle((0,0), (CR/u.kpc)*1.2, color='r', fill=False, ls='dashed')
ax.add_patch(circ)
ax.add_patch(linblad1)
ax.add_patch(linblad2)
plt.show()

duration = default_timer() - start
print 'time:'
print duration