#!/usr/bin/env python

import Orbit_Code 
reload(Orbit_Code)

orbit = Orbit_Code.Orbit_Calculator(4,25,.5,8,.4,6.5,1.8,-10.,223.) 
orbit.makeOrbit() 
#orbit.saveData("C:/Trapped_Orbital_Integrator/qp_file/")

orbit.plot(1)