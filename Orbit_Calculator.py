import astropy.units as u
import astropy.constants as const
import numpy as np
import matplotlib.pyplot as plt
from timeit import default_timer

###############################################################
# Defining some constants
###############################################################
pi = np.pi
G = const.G

###############################################################
# Defining model dependent constants
###############################################################
      
#Disk parameters (You can toggle these)
vc = 220 *u.km /u.s             # Circular velocity of the disk
RSun = 8 *u.kpc                 # Solar galactocentric radius
SigmaSun = 50 *u.Msun /u.pc**2  # Surface density at the solar radius
Rd = 2.5 *u.kpc                 # Scale length of the disk 
  
#Disk parameters (You can toggle these), Other paramaters specified in constructor
theta = 25 *u.degree


'''
The Orbit_Calculator class performs all of the calculations associated with
numerically simulating a star's orbit about the center of a spiral galaxy
and plots the results. 
'''
class Orbit_Calculator(object):

###############################################################
#Constructor
###############################################################

    def __init__(self,m1,IntTime1,CR1,epsilon1,x01,y01,vx01,vy01):
        
        #The arguments given are assigned to global variables.
        #This is obviously not the best way to use global variables in a class,
        #but in doing so the original code is more easily salvaged.
        
        global m
        global IntTime
        global CR
        global epsilon
        global x0
        global y0
        global vx0
        global vy0
        global IntTimeUnitless   
        global StepTime      
        global NSteps                  
        global T
        global OmegaCR
        
        #Instantiate global variables
        m = m1
        IntTime = IntTime1*u.Gyr
        CR = CR1*u.kpc
        epsilon = epsilon1
        x0 = x01
        y0 = y01
        vx0 = vx01
        vy0 = vy01 
        OmegaCR = vc/CR
        
        #Some time related stuff
        StepTime = 100000.*u.yr                          #Time between each calculation   
        IntTimeUnitless = (IntTime/u.yr).decompose()     #Simulation time
        NSteps = np.rint((IntTime/StepTime).decompose()) #Integer number of total steps to take
        T = np.linspace(0,IntTimeUnitless,NSteps)        #Time values
       
# Note: 
# All methods preceded by "__" are private functions that can't be called
# outside of the class 
            
###############################################################
# Method for Disk Paramaters
###############################################################

# Calculates the surface density of the disk at any radius
    def __findSurfaceDensity(self,R): 
        
        Sigma0 = SigmaSun *np.exp(RSun/Rd)
        SigmaR = Sigma0 *np.exp(-R/Rd)
        return SigmaR


###############################################################        
# Methods for spiral parameters
###############################################################

# Calculates parameter alpha
    def __findalpha(self): 
        
        global alpha
        alpha = m/np.tan(theta)
        return alpha
        
# Calculates amplitude of spiral perturbation
    def __findA(self,x,y): 
        global A
        R = np.sqrt(x**2 + y**2)
        Sigma = self.__findSurfaceDensity(R)
        A = 2. *pi *G *Sigma*epsilon*R/alpha
        return A
         
###############################################################
#Private Calculation Methods
###############################################################     

# Calculates the acceleration in this potential at coordinate x-y
    def __dvdt(self,qp,tnow):         
        
        x = qp[0] *u.kpc
        y = qp[1] *u.kpc
        A = self.__findA(x,y)
        time = tnow *u.yr
        
        # Find acceleration from logarithmic disk
        dvxD = vc**2 *x/(x**2 + y**2)
        dvyD = vc**2 *y/(x**2 + y**2)
        
        # Find acceleration from spiral
        dvSFront = -A *np.sin((m *time *vc /CR)*u.rad -m *np.arctan(y/x) -alpha *np.log(np.sqrt( x**2+y**2)/CR)*u.rad)/(x**2 +y**2)
        dvxS = dvSFront *(m*y + (-alpha-(1-np.sqrt(x**2+y**2)/Rd)*(np.tan((m *time *vc /CR)*u.rad -m *np.arctan(y/x) -alpha *np.log( np.sqrt( x**2+y**2)/CR)*u.rad)**-1))*x)
        dvyS = dvSFront *(-m*x + (-alpha-(1-np.sqrt(x**2+y**2)/Rd)*(np.tan((m *time *vc /CR)*u.rad -m *np.arctan(y/x) -alpha *np.log( np.sqrt( x**2+y**2)/CR)*u.rad)**-1))*y)
    
        # Find total acceleration
        dvxdt = ((dvxD + dvxS)/(u.km /u.s**2))
        dvydt = ((dvyD + dvyS)/(u.km /u.s**2))
        return np.array([dvxdt,dvydt])
                     
# Perform a single leapstep (t+dt), using kick-drift-kick method
    def __leapstep(self,qpnow,tnow): 
        
        dt = (StepTime/u.year)*3.15576e+07 #convert to seconds/remove units
        x = qpnow[0]
        y = qpnow[1]
        vx = qpnow[2]
        vy = qpnow[3]
        # Note: x and y are in kpc, vx and vy are in km/s
        
        a = self.__dvdt(qpnow,tnow)      # Find acceleration at this coordinate
        vx = vx -0.5 *dt *a[0]           # Advance v_x by half step
        vy = vy -0.5 *dt *a[1]           # Advance v_y by half step
        x = x +(dt*vx*3.24077928947e-17) # Advance x by full step, while converting v*dt from km to kpc
        y = y +(dt*vy*3.24077928947e-17) # Advance y by full step, while converting v*dt from km to kpc
        # x and y are in kpc, vx and vy are in km/s
        
        qpmid = np.array([x,y,vx,vy])
        a  = self.__dvdt(qpmid,tnow)     # Find a at new position and complete the velocity step
        vx = vx -0.5 *dt *a[0]           # Complete v_x step
        vy = vy -0.5 *dt *a[1]           # Complete v_y step
        
        qpnew = np.array([x,y,vx,vy,tnow]) 
        return qpnew
        
# Helper function to plot the spiral arms, arm points calcualted in polar 
# coordiantes and then converted to rectangular for plotting  
    def __plotArms(self,ax):
        
        points = 20                     #Number of points of each arm to plot
        radius = np.zeros(points)       #empty array for arm radii
        t = np.linspace(pi/24,pi/2+pi/16,points)  #array of angles 
        for i in xrange(0,m):
            radius = (CR/u.kpc)*np.exp((-m*(t)+np.pi)/alpha)
            ax.plot(radius*np.cos(t+2*np.pi*i/m),radius*np.sin(t+2*np.pi*i/m), color="purple",ls='dotted')   
            
# Convert coordinates from NR-frame to R-frame
    def __toRframe(self,qpl):  
        
        # Pull out cartesian non-rotating info
        x = qpl[:,0]
        y = qpl[:,1]
        vx = qpl[:,2]
        vy = qpl[:,3]
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
        
                                            
###############################################################
#Public Methods (Callable)
###############################################################   
    
# Calls the previously defined functions to calculate the orbit in both frames  
    def makeOrbit(self):
        
        start = default_timer()
        self.__findalpha()
    
        qp0 = np.array([x0,y0,vx0,vy0,T[0]])
        global qp
        qp = np.zeros(shape=(len(T),5))
        qp[0] = qp0
        print "Steps: %s" %str(NSteps)
        for i in range(len(T)):
            qpstep = self.__leapstep(qp0,T[i])
            qp[i] = qpstep
            qp0 = qpstep
            
        global qpR
        qpR =  self.__toRframe(qp) 
          
        duration = default_timer() - start 
        print "time: %s s" % str(duration)

# Returns qp                 
    def getqp(self):
        return qp
        
# Returns qpR        
    def getqpR(self):
        return qpR

# Saves data from non-rotating frame in dump file  
# Remember that each computer has a different file path    
    def saveData(self):
    
        filename = "qp_(m=%s)_(t=%s)_(CR=%s)_(eps=%s)_(x0=%s)_(y0=%s)_(vx0=%s)_(vy0=%s)" %(str(m),
        str(IntTime/u.Gyr),str(CR/u.kpc),str(epsilon),str(x0),str(y0),str(vx0),str(vy0))
        np.save("/Users/LBarbano/Desktop/QP_Dump/%s" % filename,qp) 
        
# Plots the orbit  
# For plot of orbit in non-rotating frame, enter 0 as the plot option
# For the plot in hte rotating frame, enter anything else        
    def plot(self,plotOption):
        
        if plotOption==0:
            qps = qp
        else:
            qps = qpR
        
        plt.close('all')         #close old plots still up
        
        fig = plt.figure(1)      #setting up the basic figure with axes and labels
        ax = fig.add_subplot(1,1,1)
        ax.set_xlabel(r'$x$ (kpc)')
        ax.set_ylabel(r'$y$ (kpc)')
        plt.axis([-13,13,-13,13])
        
        #Lindlbad Resonance stuff
        R_1o = (m+np.sqrt(2))*vc/(m*OmegaCR)
        R_1i = (m-np.sqrt(2))*vc/(m*OmegaCR)
        #for the ultraharmonic resonances
        R_2o = ((2*m)+np.sqrt(2))*vc/((2*m)*OmegaCR)
        R_2i = ((2*m)-np.sqrt(2))*vc/((2*m)*OmegaCR)

        lind1 = plt.Circle((0,0), (R_1o/u.kpc), color='g', fill=False,ls = 'dashed')
        lind2 = plt.Circle((0,0), (R_1i/u.kpc), color='g', fill=False,ls = 'dashed')
        lind1uh = plt.Circle((0,0), (R_2o/u.kpc), color='g', fill=False,ls = 'dotted')
        lind2uh = plt.Circle((0,0), (R_2i/u.kpc), color='g', fill=False,ls = 'dotted')
        ax.add_patch(lind1)
        ax.add_patch(lind2)
        ax.add_patch(lind1uh)
        ax.add_patch(lind2uh)
        
        #Plot corotation radius
        circ = plt.Circle((0,0), (CR/u.kpc), color='g', fill=False) #plotting CR radius
        ax.add_patch(circ)

        plt.plot(qps[:,0],qps[:,1], color="SlateBlue", markevery=500, marker='.', ms=8) 
        self.__plotArms(ax)
        plt.show()
        
    def doAllThings(self):
            
        self.makeOrbit()
        self.saveData()
        self.plot(1)
       
    def findEj(self):
        
        #pulling info
        x = qp[:,0]*u.kpc
        y = qp[:,1]*u.kpc
        vx = qp[:,2]*u.km/u.s
        vy = qp[:,3]*u.km/u.s
        t = qp[:,4]*u.yr
        R = np.sqrt((x**2)+(y**2))
        phi = np.arctan2(y,x)
        #finding disk potential
        disk_potential = (vc**2)*np.log(R/u.kpc)
        #finding spiral potential
        SigmaR = self.__findSurfaceDensity(CR)
        A_r = 2. *pi *G *SigmaR*epsilon*CR/alpha
        spiral_potential = A_r.to((u.km/u.s)**2)*np.cos(-alpha*np.log(R/CR)*u.rad + (m*vc*t/CR)*u.rad -m*phi)
        #calculating potential, then energy, then Ej
        potential = disk_potential + spiral_potential
        E_tot = potential + 0.5*(vx**2 + vy**2)
        L_z = R*-np.sqrt(vx**2 + vy**2)*np.sin(phi - np.arctan2(vy,vx))
        E_j = E_tot - OmegaCR*L_z
        return np.array([E_j, L_z, E_tot]) #implicit units of (km/s)**2