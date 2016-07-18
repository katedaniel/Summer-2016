#import time

import astropy.units as u
import astropy.constants as const
import numpy as np

NRun = 1 # Number of runs
NOrbit = 1 # Number of orbits per simulation

MasterOutName = "./temp_initials.txt"


################################################################################
# Defining some constants
################################################################################
pi = np.pi
G = const.G

################################################################################
# Defining model dependent constants
################################################################################
      
# Disk parameters (You can toggle these)
vc = 220. *u.km /u.s             # Circular velocity of the disk
RSun = 8. *u.kpc                 # Solar galactocentric radius
SigmaSun = 50. *u.Msun /u.pc**2  # Surface density at the solar radius
sigmaSun = 35. *u.km /u.s        # Velocity dispersion at the solar radius
Rd = 2.5 *u.kpc                  # Scale length of the disk surface density
Rs = 3.*Rd                       # Scale length for the velocity dispersion
Rp = 1. *u.kpc                   # Scale length of the disk potential
RMax = 15. *u.kpc                # Maximum radius, used to produce envelope function

'''
################################################################################
# Helpful little timer
################################################################################
def timer(start,end):
    hours, rem = divmod(end-start, 3600)
    minutes, seconds = divmod(rem, 60)
    return "{:0>2}:{:0>2}:{:05.2f}".format(int(hours),int(minutes),seconds)
'''

################################################################################
# Definitions to calculate stellar paremeters
################################################################################

# Calculates the surface density of the disk at any radius
def findVelocityDispersion(R):         
    sigma0 = sigmaSun *np.exp(RSun/Rs)
    sigmaR = sigma0 *np.exp(-R/Rs)
    return sigmaR

# Calculates the surface density of the disk at any radius
def findSurfaceDensity(R):     
    Sigma0 = SigmaSun *np.exp(RSun/Rd)
    SigmaR = Sigma0 *np.exp(-R/Rd)
    return SigmaR

# Find radius associated with a circular orbit for star with energy E (km/s)^2    
def findRE(E):
    RE = Rp *np.exp((E - vc**2/2.)/vc**2)
    return RE

# Find guiding center radius R_g    
def findRL(Lz):
    RL = Lz/vc
    return RL
    
# Find energy for star in axisymmetric potential
def findE(R,v):
    disk_potential = (vc**2)*np.log(R/Rp)
    kineticE = 0.5*v**2
    E = disk_potential + kineticE
    return E

# Find the random energy of a star -- returns (km/s)^2
def findEran(E,Lz):
    RL = findRL(Lz) # Find guiding center radius
    E_tot_c = findE(RL,vc)
    ## Calculating random energy [E_ran = E_tot - E_c(R_L)] *See expression below eqn.5 in DW15
    E_ran = E - E_tot_c
    return E_ran

################################################################################
# Monte Carlo - Rejection Method
################################################################################

# Calculate the distribution function f(E,Lz)
def findf(E,Lz):
    RE = findRE(E)
    sigmaR = findVelocityDispersion(RE)
    Sigma = findSurfaceDensity(RE)
    Lc = RE *vc
    Omega = vc/RE
    FrontStuff = Sigma/(np.sqrt(2)*pi*sigmaR**2)
    eStuff = Omega*(Lz-Lc)/sigmaR**2
    f = FrontStuff *np.exp(eStuff)
    return f

# Find envelope function (fMax)
# The sweet spot (SS) for maximizing the normalization is at R = 3 Rd
# Maximize the envelope function at SS
def findfMax(E,Lz):
    RE = findRE(E)
    sigmaR = findVelocityDispersion(RE) 
    Sigma = findSurfaceDensity(RE)
    Lc = RE *vc
    Omega = vc/RE
    FrontStuff = Sigma/(np.sqrt(2)*pi*sigmaR**2)
    eStuff = Omega*(Lz-Lc)/sigmaR**2
    fMax = FrontStuff #*np.exp(eStuff)
    return fMax

### Define max values for even distribution of E and Lz
EMax = vc**2*np.log(RMax/Rp) + 0.5*(vc)**2
LzMax = RMax*vc
M = 1.

def getMCqp0(): # Define initial conditions evenly distributed within f_New    
    nOK = 0 # An approved set of inital conditions has not yet been produced
    while (nOK < 1):
        # Initial random values
        iLz = LzMax *np.random.random()
        iE = EMax *np.random.random()
        # Declare the value of the DF f(x)
        fx = findf(iE,iLz)
        # Find enveloope function
        fMax = findfMax(EMax,LzMax)
        gx = fMax*M
        # Determine acceptance probability
        AcceptProb = fx/gx
        # Provide a uniform random variable
        utest = np.random.random()
        # If utest < AcceptProb then accept these values for E and Lz
        if utest < AcceptProb:
            RL = findRL(iLz)/u.kpc
            if (RL > 4) and (RL < 15):
                Eran = findEran(iE,iLz) # Find random energy(km/s)^2
                if (Eran/(u.km/u.s)**2 > 0):
                    vran = np.sqrt(2.*Eran) # Find amplitude of random velocity (km/s)
                    if (vran < 2.*findVelocityDispersion(3.*Rd)):
                        ranglec = 2.*pi *np.random.random() # Produce direction for position vector
                        xc = RL *np.cos(ranglec)
                        yc = RL *np.sin(ranglec)
                        # Define coords related to elliptical epicycle
                        ae = 1./np.sqrt(2.)
                        be = 0.5
                        rangleran = 2.*pi *np.random.random() # Direction for position in epicycle
                        xe = ae *np.cos(rangleran)
                        ye = ae *np.sin(rangleran)
                        Rran = np.sqrt(xe**2 + ye**2)
                        vanglerane = np.arctan2(-be**2 *xe, ae**2 *ye) # slope tangent to ellips in xe-ye
                        vangleran = vanglerane -(ranglec -pi/2.) # slope tangent to ellipse in x-y
                        xran = Rran *np.cos(rangleran -ranglec +pi/2.)
                        yran = Rran *np.sin(rangleran -ranglec +pi/2.)
                        x0 = xc + xran
                        y0 = yc + yran
                        # Figure out the velocity
                        vcx = vc *np.cos(ranglec +pi/2.)
                        vcy = vc *np.sin(ranglec +pi/2.)            
                        vranx = vran *np.cos(vangleran) # x-component of random velocity (km/s)
                        vrany = vran *np.sin(vangleran) # y-component of random velocity (km/s)
                        vx0 = vcx + vranx # initial x-component of velocity vector in N-frame (km/s)
                        vy0 = vcy + vrany # initial y-component of velocity vector in N-frame (km/s)            
                        R = np.sqrt(x0**2 +y0**2)
                        if  (R > 2) and (R < 15):
                            x0 = x0
                            y0 = y0
                            # Get rid of explicit units
                            vx0 = vx0/(u.km/u.s)
                            vy0 = vy0/(u.km/u.s)
                            qp0 = np.array([x0,y0,vx0,vy0])
                            nOK = 1
    return qp0    

nRun = 1
#starttime = time.time()
'''
while nRun < NRun+1:
    AOut = []
    nOrbit = 1
    while nOrbit < NOrbit+1:
        qp0 = getMCqp0()
        AOut.append([qp0[0],qp0[1],qp0[2],qp0[3],nRun,nOrbit])
        rangle = np.arctan2(qp0[1],qp0[0])
        vangle = np.arctan2(qp0[3],qp0[2])
        alph = vangle - rangle
        vtot = np.sqrt(qp0[2]**2 + qp0[3]**2)
        vr = np.round(vtot*np.cos(alph),1)
        vphi = np.round(vtot*np.sin(alph) - 220.,1)
        R = np.round(np.sqrt(qp0[0]**2 + qp0[1]**2),1)
        vran = np.round(np.sqrt(vr**2 +vphi**2),1)
        print 'R = ',R, '| v_R = ', vr, '| v_phi = ', vphi, '| v_ran = ', vran
        #print "Orbit # ", nOrbit, "| Run # ", nRun, "| Time Elapsed: ", timer(starttime,time.time()), "| ", round(100.*(float(nRun-1)*float(NOrbit) + (nOrbit-1))/float(NRun*NOrbit),1), "% Finished |"
        nOrbit = nOrbit+1
    nRun = nRun +1
'''
AOut = []
qp0 = getMCqp0()
AOut.append([qp0[0],qp0[1],qp0[2],qp0[3]])
initials = AOut[0]
np.savetxt(MasterOutName,initials, delimiter="", fmt="%s", newline=" ")
