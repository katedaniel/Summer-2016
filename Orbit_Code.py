'''
Authors: Luke Barbano and Noah Lifset

Description: 
    This file implements the Orbit_Calculator class, which performs all of the 
calculations associated with numerically simulating a star's orbit about the 
center of a spiral galaxy and plots the results. 
    The class takes in various parameters relating to the spiral as well as the 
star's initial conditions. It can then calculate the positions of the star's 
orbit for some time interval, plot the results in the non-rotating frame, 
translate the positions to the rotating frame of the spiral, plot the orbit in 
the rotating frame, and save the position/velocity information from the non-
rotating frame to a dump file in a specified file path. The accompanying file 
"OrbitDemonstration.py" provides a walkthrough of how to use this class.  
    As an aside, this class is not super well written by computer science standards
given its reliance on global variables in lieu of assigning attributes to the 
instance of the class; however, this class implementation best salvaged the 
previous version of the code and might be slightly faster since the self instance 
of the class doesn't have to constanly access the many, frequently used variables. 

Note: If any changes are made to the Orbit_Calculator class, make sure to restart the 
kernel so the changes will be implemented when using the class.

Class Summary:
    
--------------------------------------------------------------------------------
class Orbit_Calculator(__builtin__.object)
       
Constructor:     
     __init__(self, m, IntTime, CR, epsilon, x0, y0, vx0, vy0)
        m = number of spiral arms
        IntTime = Duration of simulation (units implicit gigayears)
        CR = Corotation radius (units implicit kiloparsecs)
        epsilon = epsilon of spiral
        x0,y0,vx0,vy0 = intial x,y,vx,and vy of star (implicit units kpc and km/s respectively)
    
Public (Callable) Methods:
    
    --makeOrbit(self)
        -calculates the orbit
        -creates global numpy arrays for qp and qpR

    --getqp(self)
        -returns numpy array qp
        
    --getqpR(self)
        -returns numpy array qpR        
                       
    --saveData(self)  
        -save qp to a file in a designated filepath 
        
    --plot(self, plotOption)
        -plots the orbit
        -if plotOption == 0, plots in the non-rotating frame
        -if plotOption is anything but 0, plots in the rotating frame
         
    --doAllThings(self)
        -calls makeOrbit(), saveData(), and plot(1)
        -calculates the orbit, saves qp, and plot orbit in rotating frame
    
    --findEj(self)
        -calculates Ej (Jacobi Integral)
        -returns numpy array Ej
        
    --Capture(self)
        -calculates capture criteria for guiding center
        -returns numpy array of lambdas (capture criteria) 
    
_______________________________________________________________________________    
'''
import astropy.units as u
import astropy.constants as const
import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt
from timeit import default_timer
from numpy import arange
from numpy import meshgrid
from mpl_toolkits.mplot3d import Axes3D

################################################################################
# Defining some constants
################################################################################
pi = np.pi
G = const.G

################################################################################
# Defining model dependent constants
################################################################################
      
#Disk parameters (You can toggle these)
vc = 220 *u.km /u.s             # Circular velocity of the disk
RSun = 8 *u.kpc                 # Solar galactocentric radius
SigmaSun = 50 *u.Msun /u.pc**2  # Surface density at the solar radius
Rd = 2.5 *u.kpc                 # Scale length of the disk 
  

class Orbit_Calculator(object):

################################################################################
#Constructor
################################################################################

    def __init__(self,m1,theta1,IntTime1,CR1,epsilon1,x01,y01,vx01,vy01):
        
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
        global theta 
        
        #Assign global variables
        m = int(m1)
        theta = theta1*u.degree
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
        
        self.__findalpha()
        self.__findhcr()
# Note: 
# All methods preceded by "__" are private functions that can't be called
# outside of the class 
            
################################################################################
# Method for Disk Parameters
################################################################################

# Calculates the surface density of the disk at any radius
    def __findSurfaceDensity(self,R): 
        
        Sigma0 = SigmaSun *np.exp(RSun/Rd)
        SigmaR = Sigma0 *np.exp(-R/Rd)
        return SigmaR


################################################################################        
# Methods for spiral parameters
################################################################################

# Calculates parameter alpha
    def __findalpha(self): 
        
        global alpha
        alpha = m/np.tan(theta)
        
# Calculates amplitude of spiral perturbation
    def __findA(self,R): 
        Sigma = self.__findSurfaceDensity(R)
        A = 2. *pi *G *Sigma*epsilon*R/alpha
        return A
        
# Calculates hcr, Jacobi integral for a star at corotation with no spirals
    def __findhcr(self):
        disk_potential = (vc**2)*np.log(CR/u.kpc)
        E_tot = disk_potential + 0.5*(vc**2)
        L_z = CR*vc
        global hcr
        hcr = (E_tot - OmegaCR*L_z) #units of (km/s)^2     
        
                
################################################################################
#Private Calculation Methods
################################################################################    

# Calculates the acceleration in this potential at coordinate x-y
    def __dvdt(self,qp,tnow):         
        
        x = qp[0] *u.kpc
        y = qp[1] *u.kpc
        time = tnow *u.yr
        R = (x**2 + y**2)
        A = self.__findA(np.sqrt(R))
        
        # Find acceleration from logarithmic disk
        dvxD = vc**2 *x/R
        dvyD = vc**2 *y/R
        
        # Find acceleration from spiral
        var1 = (m *time *vc /CR)*u.rad -m *np.arctan(y/x) -alpha *np.log(np.sqrt(R)/CR)*u.rad
        var2 = (-alpha-(1-np.sqrt(R)/Rd)*(np.tan(var1)**-1))
        dvSFront = -A *np.sin(var1)/(R)
        dvxS = dvSFront *(m*y + var2*x)
        dvyS = dvSFront *(-m*x + var2*y)
    
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
        
        return np.array([x,y,vx,vy,tnow]) 
        
# Helper function to plot the spiral arms, arm points calcualted in polar 
# coordinates and then converted to rectangular for plotting  
    def __plotArms(self,ax):
        
        points = 20                     #Number of points of each arm to plot
        radius = np.zeros(points)       #empty array for arm radii
        t = np.linspace(pi/24,pi/2+pi/16,points)  #array of angles 
        for i in xrange(0,m):
            radius = (CR/u.kpc)*np.exp((-m*(t)+np.pi)/alpha)
            ax.plot(radius*np.cos(t+2*np.pi*i/m),radius*np.sin(t+2*np.pi*i/m), color="purple",ls='dotted')   
        return
          
# Convert coordinates from NR-frame to R-frame
    def __toRframe(self,qpl):  
        
        # Pull out cartesian non-rotating info
        x = qpl[:,0]
        y = qpl[:,1]
        vx = qpl[:,2]
        vy = qpl[:,3]
        t = qpl[:,4]*u.yr
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
        
                                            
################################################################################
#Public Methods (Callable)
################################################################################  
    
# Calls the previously defined functions to calculate the orbit in both frames  
    def makeOrbit(self):
        
        start = default_timer()
    
        qp0 = np.array([x0,y0,vx0,vy0,T[0]])
        global qp
        qp = np.zeros(shape=(len(T),5))
        qp[0] = qp0
        print "Steps: %s" %str(NSteps)
        for i in range(len(T)):
            qpstep = qp0 = self.__leapstep(qp0,T[i])
            qp[i] = qpstep
            
        global qpR
        qpR =  self.__toRframe(qp) 
          
        duration = default_timer() - start 
        print "time: %s s" % str(duration)
        return

# Returns qp                 
    def getqp(self):
        return qp
        
# Returns qpR        
    def getqpR(self):
        return qpR

# Returns NSteps                 
    def getNSteps(self):
        return NSteps
        
# Sets qp                 
    def setqp(self,qps):
        global qp
        global qpR
        qp = qps
        qpR = self.__toRframe(qp)
        return
        
# Saves data from non-rotating frame in dump file  
# Remember that each computer has a different file path    
    def saveData(self,filepath,):
        filename = "qp_(m=%s)_(th=%s)_(t=%s)_(CR=%s)_(eps=%s)_(x0=%s)_(y0=%s)_(vx0=%s)_(vy0=%s)" %(str(m),
        str(theta/u.degree),str(IntTime/u.Gyr),str(CR/u.kpc),str(epsilon),str(x0),str(y0),str(vx0),str(vy0))
        np.save(filepath + filename,qp) 
        return
        
# Plots the orbit  
# For plot of orbit in non-rotating frame, enter 0 as the plot option
# For the plot in the rotating frame, enter anything else        
    def plot(self,plotOption):
        
        plt.close('all')         #close old plots still up
        
        fig = plt.figure(1)      #setting up the basic figure with axes and labels
        ax = fig.add_subplot(1,1,1)
        ax.set_xlabel(r'$x$ (kpc)')
        ax.set_ylabel(r'$y$ (kpc)')
        plt.axis([-13,13,-13,13])
        ax.set_aspect('equal', 'datalim')
        
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
        
        self.__plotArms(ax)
        plt.show()
        
        #Plot corotation radius
        circ = plt.Circle((0,0), (CR/u.kpc), color='g', fill=False) #plotting CR radius
        ax.add_patch(circ)
        
        #Plot the capture region
        if plotOption != 0:
            delta = 0.025
            x_range = arange(-13.0, 13.0, delta)
            y_range = arange(-13.0, 13.0, delta)
            X, Y = meshgrid(x_range,y_range)
            R = np.sqrt((X**2) + (Y**2))
            phi = np.arctan2(Y,X)
            #find Phi_eff_min
            A_CR = self.__findA(CR).to((u.km/u.s)**2)
            phi_min = (hcr - A_CR)/((u.km/u.s)**2)
            phi_max = (hcr + A_CR)/((u.km/u.s)**2)
            #defining the contour equation
            A = self.__findA(R*u.kpc).to((u.km/u.s)**2)
            spiral_potential = A_CR*np.cos(-alpha*np.log(R*u.kpc/CR)*u.rad -m*phi*u.rad)
            disk_potential = (vc**2)*np.log(R)
            potential = spiral_potential + disk_potential
            func = (0.5*(OmegaCR**2)*(CR**2) - (OmegaCR**2)*CR*(R*u.kpc) + potential)/((u.km/u.s)**2)
            #plotting contour
            plt.contourf(X,Y,func,[phi_min,phi_max],colors='gray',alpha=0.3)
            
        if plotOption==0:
            qps = qp
        elif plotOption==1:
            qps = qpR
            [Rgx,Rgy] = self.findRg()
            plt.plot(Rgx,Rgy,color="black", markevery=500, marker='.', ms=8)
        else:
            return fig, ax  
        plt.plot(qps[:,0],qps[:,1], color="SlateBlue", markevery=500, marker='.', ms=8) 
        
        return fig, ax
    
        
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
        A = self.__findA(R).to((u.km/u.s)**2)
        spiral_potential = A*np.cos(-alpha*np.log(R/CR)*u.rad + (m*OmegaCR*t)*u.rad -m*phi)
        #calculating potential, then energy, then Ej
        potential = disk_potential + spiral_potential
        E_tot = potential + 0.5*(vx**2 + vy**2)
        L_z = R*-np.sqrt(vx**2 + vy**2)*np.sin(phi - np.arctan2(vy,vx))
        E_j = E_tot - OmegaCR*L_z
        return np.array([E_j, L_z, E_tot]) #implicit units of (km/s)**2
        
        
    def Capture(self):
        
        #pulling info out of qp
        x = qp[:,0]*u.kpc
        y = qp[:,1]*u.kpc
        vx = qp[:,2]*u.km/u.s
        vy = qp[:,3]*u.km/u.s
        R = np.sqrt(x**2 + y**2)
        phi = np.arctan2(y,x)
        #pulling some constants
        Ej__ = self.findEj()
        Ej_ = Ej__[0]
        Ej = Ej_[0]
        A_CR = self.__findA(CR).to((u.km/u.s)**2)
        #finding Rg
        R_g = R*-np.sqrt(vx**2 + vy**2)*np.sin(phi - np.arctan2(vy,vx))/vc
        #finding E_random
        E_ran = 0.5*((np.sqrt(vx**2 + vy**2)*np.cos(phi - np.arctan2(vy,vx)))**2 + (-np.sqrt(vx**2 + vy**2)*np.sin(phi - np.arctan2(vy,vx))-vc)**2)
        #finding Lambda_c
        Lam_c = (Ej - hcr/((u.km/u.s)**2))/(A_CR/((u.km/u.s)**2))
        #finding Lambda_nc,2
        Lam_nc2 = Lam_c - ((R_g/CR)*(E_ran/A_CR))
        return np.array(Lam_nc2)
        
    def Phi_eff(self):
        
        #pulling info out of qp
        x = qp[:,0]*u.kpc
        y = qp[:,1]*u.kpc
        vx = qp[:,2]*u.km/u.s
        vy = qp[:,3]*u.km/u.s
        t = qp[:,4]*u.yr
        R = np.sqrt(x**2 + y**2)
        phi = np.arctan2(y,x)
        A = self.__findA(R).to((u.km/u.s)**2)
        #finding potential
        disk_potential = (vc**2)*np.log(R/u.kpc)
        spiral_potential = A*np.cos(-alpha*np.log(R/CR)*u.rad + m*phi)
        potential = disk_potential + spiral_potential
        #put it together
        phi_eff = potential - 0.5*(OmegaCR*R)**2
        return phi_eff
        
    def test(self):
        
        delta = 0.25
        x_range = arange(-13.0, 13.0, delta)
        y_range = arange(-13.0, 13.0, delta)
        X, Y = meshgrid(x_range,y_range)
        R = np.sqrt((X**2) + (Y**2))*u.kpc
        phi = np.arctan2(Y,X)
        A = self.__findA(R).to((u.km/u.s)**2)
        disk_potential = (vc**2)*np.log(R/u.kpc)
        spiral_potential = A*np.cos(-alpha*np.log(R/CR)*u.rad - m*phi*u.rad)
        potential = disk_potential + spiral_potential
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.plot_wireframe(X,Y,potential - 0.5*(OmegaCR*R)**2)#plots phi_eff
        return

# Calculates position of guiding radius in rotating frame        
    def findRg(self):
        
        #pulling info out of qp
        x = qp[:,0]*u.kpc
        y = qp[:,1]*u.kpc
        vx = qp[:,2]*u.km/u.s
        vy = qp[:,3]*u.km/u.s
        t = qp[:,4] *u.yr
        R_g = (np.sqrt(x**2 + y**2)*-np.sqrt(vx**2 + vy**2)*np.sin(np.arctan2(y,x) - np.arctan2(vy,vx))/vc)
        phi = np.arctan2(y,x)
        phiR= phi - (t*OmegaCR).decompose() *u.rad
        xR = R_g *np.cos(phiR)
        yR = R_g *np.sin(phiR)
        return np.array([xR,yR])
        
    def Poincare(self):
        plt.close('all')
        yspline = interpolate.splrep(qp[:,4], qpR[:,1], s=0)
        roots = interpolate.sproot(yspline)
        if len(roots)==0:
            return 
        xspline = interpolate.splrep(qp[:,4], qpR[:,0], s=0)
        vxspline = interpolate.splrep(qp[:,4], qpR[:,2], s=0)
        x = interpolate.splev(roots, xspline)
        vx = interpolate.splev(roots, vxspline)
        plt.scatter(x,vx)
        plt.show()
        return 
        
        
       
        