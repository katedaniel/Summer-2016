import Orbit_Code 
reload(Orbit_Code)

orbit = Orbit_Code.Orbit_Calculator(4,25,.1,8,.4,7.8,0.,3.,230.) 
orbit.makeOrbit() 
orbit.saveData("/Users/LBarbano/Desktop/QP_Dump/")