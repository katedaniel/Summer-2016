import Orbit_Code
reload(Orbit_Code)

orbit = Orbit_Code.Orbit_Calculator(4,1.5,7.,0.3,6.8,0.,3.,223.)

orbit.makeOrbit()

orbit.plot(1)

orbit.saveData()