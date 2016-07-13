import Orbit_Code 
reload(Orbit_Code)

orbit = Orbit_Code.Orbit_Calculator(4,25,.5,8,.3,7.5,2.,30.,230.) 
orbit.makeOrbit() 
#orbit.saveData("C:/Users/Noah/Documents/GitHub/Trapped_Orbital_Integrator/qp_file/")
#tck = orbit.Poincare()