# Simple optimization example - minimize constrained Betts function

import scipy.optimize as sopt
import numpy as np
import rref


w = np.concatenate((np.ones((1, 33)), 0.01*np.ones((1, 8)), 10*np.ones((1, 9))), axis=1)
P = 2*np.diag(w)

lb = np.concatenate((np.zeros((1, 33)), -6*np.ones((1, 17))), axis=1)
ub = np.concatenate((np.ones((1, 33)), 6*np.ones((1, 17))), axis=1)
bounds = np.concatenate((np.transpose(lb), np.transpose(ub)), axis=1)
bounds_tuple = tuple(map(tuple, bounds))

# linear equality constraint of the form A_eq x = B_eq
A_eq = np.load('A_eq.npy')
B_eq = np.load('B_eq.npy')
x0 = np.load('x0.npy')

opts = {'maxiter': 10000, 'ftol': 1E-3, 'disp': True, 'eps': 1E-6}

# reduce matrix A_eq
A_eq_reg, index_pivots, index_rearranged = rref.rref(A_eq)
rank = len(index_pivots)
A_eq_reg = A_eq_reg[0:rank, :]
B_eq_reg = B_eq[index_rearranged]
B_eq_reg = B_eq_reg[0:rank]

# set up the constraints
lc_eq = sopt.LinearConstraint(A_eq_reg, lb=B_eq_reg, ub=B_eq_reg)
cons = [lc_eq]

result = sopt.minimize(obj_funct, x0, method='SLSQP', constraints=cons, bounds=bounds_tuple, options=opts)
x_opt = result.x
