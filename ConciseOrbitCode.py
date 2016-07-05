#!/usr/bin/env python

import Orbit_Code 
reload(Orbit_Code)

orbit = Orbit_Code.Orbit_Calculator(4,25,.1,8,.4,7.8,0.,-10.,223.) 
orbit.makeOrbit() 

#orbit.saveData("C:/Trapped_Orbital_Integrator/qp_file/")
#orbit.plot(1)
#tck = orbit.Poincare()