import Orbit_Code
reload(Orbit_Code)

from Orbit_Calculator import Orbit_Calculator

orbit = Orbit_Code.Orbit_Calculator(4,0.1,8.,0.3,7.8,0.,3.,223.)

orbit.makeOrbit()