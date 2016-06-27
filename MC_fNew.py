import astropy.units as u
import astropy.constants as const
import numpy as np

RMax = 15. *u.kpc # Maximum radius, used to produce envelope function
N = 2 # Number of runs
n = 10 # Number of orbits per simulation

################################################################################
# Defining some constants
################################################################################
pi = np.pi
G = const.G

################################################################################
# Defining model dependent constants
################################################################################
      
# Disk parameters (You can toggle these)
vc = 220 *u.km /u.s             # Circular velocity of the disk
RSun = 8 *u.kpc                 # Solar galactocentric radius
SigmaSun = 50 *u.Msun /u.pc**2  # Surface density at the solar radius
sigmaSun = 35 *u.km /u.s        # Velocity dispersion at the solar radius
Rd = 2.5 *u.kpc                 # Scale length of the disk surface density
Rs = 3.*Rd                      # Scale length for the velocity dispersion
Rp = 1. *u.kpc                  # Scale length of the disk potential


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
    print 'E = ', E
    E_tot_c = findE(RL,vc)
    print 'E_tot_c = ', E_tot_c
    ## Calculating random energy [E_ran = E_tot - E_c(R_L)] *See expression below eqn.5 in DW15
    E_ran = E - E_tot_c
    print 'E_ran = ', E_ran
    return E_ran

# Calculate the distribution function f(E,Lz)
def findf(E,Lz):
    RE = findRE(E)
    sigmaR = findVelocityDispersion(RE)
    Sigma = findSurfaceDensity(RE)
    Lc = RE *vc
    Omega = vc/RE
    FrontStuff = Sigma/(np.sqrt(2)*pi*sigmaR**2)
    eStuff = Omega*(Lc-Lz)/sigmaR**2
    f = FrontStuff *np.exp(eStuff)
    return f

################################################################################
# Monte Carlo - Rejection Method
################################################################################

### Define max values for even distribution of E and Lz
EMax = vc**2*np.log(RMax/Rp) + 0.5*(vc+findVelocityDispersion(RMax))**2
LzMax = RMax*vc

### Find envelope function
# The sweet spot (SS) for maximizing the normalization is at R = 3 Rd
# Maximize the envelope function at SS
sigmaRSS = findVelocityDispersion(3.*Rd) 
SigmaSS = findSurfaceDensity(3.*Rd)
FrontStuff = SigmaSS/(np.sqrt(2)*pi*sigmaRSS**2)
# Set the exponential part to 1
eStuff = np.log(1)
fMax = FrontStuff *np.exp(eStuff)
M = 1.
gx = fMax*M

def getMCqp0(): # Define initial conditions evenly distributed within f_New    
    nOK = 0 # An approved set of inital conditions has not yet been produced
    while (nOK < 1):
        # Initial random values
        iLz = LzMax *np.random.random()
        iE = EMax *np.random.random()
        # Declare the value of the DF f(x)
        fx = findf(iE,iLz)
        # Determine acceptance probability
        AcceptProb = fx/gx
        # Provide a uniform random variable
        utest = np.random.random()
        # If u < AcceptProb then accept these values for E and Lz
        if utest < AcceptProb:
            rangle = 2.*pi *np.random.random() # Produce direction for position vector
            vangle = 2.*pi *np.random.random() # Produce direction for velocity vector
            Eran = findEran(iE,iLz) # Find random energy(km/s)^2
            vran = np.sqrt(2.*Eran) # Find amplitude of random velocity (km/s)
            print '|v_ran| = ', vran
            vranx = vran *np.cos(vangle) # x-component of random velocity (km/s)
            vrany = vran *np.sin(vangle) # y-component of random velocity (km/s)
            vcx = vc *np.cos(rangle) # x-component of circular velocity (km/s)
            vcy = vc *np.sin(rangle) # y-component of circualar velocity (km/s)
            vx0 = vcx + vranx # initial x-component of velocity vector in N-frame (km/s)
            vy0 = vcy + vrany # initial y-component of velocity vector in N-frame (km/s)
            alph = rangle - vangle # angle between position and velocity vectors
            v_tot = np.sqrt(vx0**2 + vy0**2) # amplitude of velocity vector in N-frame (km/s)
            vphi = -v_tot*np.sin(alph) # azimuthal velocity
            R = iLz/vphi
            x0 = R *np.cos(rangle)/(u.kpc)
            y0 = R *np.sin(rangle)/(u.kpc)
            # Get rid of explicit units
            vx0 = vx0/(u.km/u.s)
            vy0 = vy0/(u.km/u.s)
            qp0 = [x0,y0,vx0,vy0]
            nOK = 1
    return qp0    

qp0= getMCqp0()
#print qp0