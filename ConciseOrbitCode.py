import Orbit_Code 
reload(Orbit_Code)

orbit = Orbit_Code.Orbit_Calculator(4,25,.01,8,.4,7.8,0.,3.,230.) 
orbit.makeOrbit() 
orbit.saveData("C:/Users/Noah/Documents/GitHub/Trapped_Orbital_Integrator/qp_file/")
tck = orbit.Poincare()