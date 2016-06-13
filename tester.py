import matplotlib.pyplot
from numpy import arange
from numpy import meshgrid
import numpy as np

delta = 0.025
xrange = arange(-13.0, 13.0, delta)
yrange = arange(-13.0, 13.0, delta)
X, Y = meshgrid(xrange,yrange)
R = np.sqrt(X**2 + Y**2)
phi = np.arctan2(Y,X)


matplotlib.pyplot.contourf(X, Y, 4. - R*phi - R**2, [1.,8],colors='gray',alpha=0.5)
matplotlib.pyplot.show()