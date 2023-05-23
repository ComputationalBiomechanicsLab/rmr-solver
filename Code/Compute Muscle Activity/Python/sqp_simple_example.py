# Test script for solving easy optimization problem with FMINCON and SQP
# Credits: https://www.youtube.com/watch?v=8uDpy-pbJr0
import numpy as np
import scipy.optimize as sopt
import time


# define function to calculate volume of box
def calcVolume(x):
    length = x[0]
    width = x[1]
    height = x[2]
    volume = length*width*height
    return volume


# define function to calculate volume of box
def calcSurface(x):
    length = x[0]
    width = x[1]
    height = x[2]
    surfaceArea = 2*(length*width + length*height + height*width)
    return surfaceArea


# define objective function for optimization
def obj_fnct(x, args):
    return -calcVolume(x)


# set initial guess values for box dimensions
length0 = 1
width0 = 1
height0 = 1

x0 = np.array([length0, width0, height0])

# set up the constraints
nlc = sopt.NonlinearConstraint(calcSurface, 0, 10)

# call solver to minimize the objective function given constraints
start = time.time()
aux = 1.4
result = sopt.minimize(obj_fnct, x0=x0, args=aux, method='SLSQP', constraints=nlc)
end = time.time() - start
print(end)
xopt = result.x

volumeOpt = calcVolume(xopt)

