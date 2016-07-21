#!/usr/bin/env python

'''
This file demonstrates how to utilize the Orbit_Calculator Class.
Start by importing the Orbit_Calculator class from the file Orbit_Calculator.
Once it has been imported, type help(Orbit_Calculator) to get a summary of the 
class's attributes and methods. See the command line for the class summary 
provided by the help() function. Note that this help step is not necessary, merely
helpful when using this class at the beginning.
'''
import Orbit_Code 
reload(Orbit_Code)

#help(Orbit_Code.Orbit_Calculator)
'''
Instantiate an Orbit_Calculator object with the appropriate parameters. The 
instance of our Orbit_Calculator will be called "orbit"
The arguments are respectively:
1: number of spiral arms
2: Duration of simulation (units implicit gigayears)
3: Corotation radius (units implicit kiloparsecs)
4: epsilon of spiral
5,6,7,8: x0,y0,vx0,vy0 (implicit units of kpc and km/s respectively)
'''
orbit = Orbit_Code.Orbit_Calculator(4,25,.1,8,.4,7.8,0.,3.,230.) 

'''
Calculate the orbit. Once you've calculated the orbit, there is no need to 
recalculate it unless you want to run a simulation with different parameters.
Note that makeOrbit() must 
'''
orbit.makeOrbit() 

'''
Plot the orbit in the non-rotating or rotating frame. The plot method takes an 
argument of 0 to plot the star's orbit in the non-rotating frame and anything else for 
the rotating frame.
'''
fig,ax = orbit.plot(1) 

'''
To plot the orbit in the rotating frame, it is not necessary to redo the 
calculations. Merely use the above method with an argument of 1 (or anything 
else i.e. 2.0, 100, 'hey' etc). The orbit object will close previous plots and 
replot the plot in the rotating frame. It makes the most sense to plot subsequent 
plots from the command line, otherwise only the last plot command would show.
'''
#orbit.plot(0)

'''
We can also access qp and qpR with the following methods. qp and qpR will be 
stored in a and b respectively
'''  

#a = orbit.getqp()
#b = orbit.getqpR()

'''
Now we can save the data in the non-rotating frame to our dump folder by using the 
saveData method. Specifcy a filepath as an argument when saving the data
'''
filename = "/Users/LBarbano/Desktop/QP_Dump/"
orbit.saveData(filename)

'''
Finally, we can do all of these things simply with the method doAllThings()
This method will call makeOrbit(), plot(1), and saveData() in succession. 
This is helpful when we just want to run the simulation, plot the orbit in the 
rotating frame, and save qp. Type the following code into the command line to
try it out. It should do the same thing as the previous code, except that it 
plots the orbit in the rotating frame and doesn't store qp and qpR in a and b. 
Not that you must have already run a file that imported Orbit_Calculator in order
to use Orbit_Calculator from the command line.
'''
#orbit1 = Orbit_Calculator(4,.01,8,.4,7.8,0.,3.,230.)
#orbit1.doAllThings()
