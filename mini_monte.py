import numpy as np

N = 1000000. #number of points
array = np.random.rand(N,2)

circle_counter = 0

for point in array:
    R = np.sqrt((point[0]**2) + (point[1]**2))
    if R <= 1.:
        circle_counter += 1

pi = (circle_counter/N)*4.

print pi