import Orbit_Code
reload(Orbit_Code)

orbit = Orbit_Code.Orbit_Calculator(4,25.,.1,8.,0.3,7.8,0.,3.,223.)

orbit.makeOrbit()


print orbit.getqp()
