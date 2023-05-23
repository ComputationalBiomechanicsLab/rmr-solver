# Simple optimization example - minimize constrained Betts function

import scipy.optimize as sopt
import numpy as np
import rref
import osqp
from scipy import sparse


def obj_funct(x):
    w = np.concatenate((np.ones((1, 33)), 0.01*np.ones((1, 8)), 10*np.ones((1, 9))), axis=1)
    cost = w.dot(np.square(x))
    return cost


lb = np.concatenate((np.zeros((1, 33)), -6*np.ones((1, 17))), axis=1)
ub = np.concatenate((np.ones((1, 33)), 6*np.ones((1, 17))), axis=1)
bounds = np.concatenate((np.transpose(lb), np.transpose(ub)), axis=1)
bounds_tuple = tuple(map(tuple, bounds))

# linear equality constraint of the form A_eq x = B_eq
A_eq = np.load('A_eq.npy')
B_eq = np.load('B_eq.npy')
x0 = np.load('x0.npy')

opts = {'maxiter': 10000, 'ftol': 1E-3, 'disp': True, 'eps': 1E-6}

# # reduce matrix A_eq
# A_eq_reg, index_pivots, index_rearranged = rref.rref(A_eq)
# rank = len(index_pivots)
# A_eq_reg = A_eq_reg[0:rank, :]
# B_eq_reg = B_eq[index_rearranged]
# B_eq_reg = B_eq_reg[0:rank]

# A_eq = np.delete(A_eq, [6, 7, 8], 0)
# B_eq = np.delete(B_eq, [6, 7, 8], 0)

## SOLVE with scipy.optimize.minimize
# set up the constraints
lc_eq = sopt.LinearConstraint(A_eq, lb=B_eq*0.99, ub=B_eq*1.01)
cons = [lc_eq]

result = sopt.minimize(obj_funct, x0, method='SLSQP', constraints=cons, bounds=bounds_tuple, options=opts)
x_opt = result.x

## SOLVE with OSQP
w = np.concatenate((np.ones((1, 33)), 0.01*np.ones((1, 8)), 10*np.ones((1, 9))), axis=1)
P = sparse.csc_matrix(2*np.diagflat(w))
A_eq_osqp = sparse.csc_matrix(np.concatenate((A_eq, np.eye(50)), axis=0))
lb_osqp = np.concatenate((np.transpose(np.atleast_2d(B_eq)), np.transpose(lb)), axis=0)
ub_osqp = np.concatenate((np.transpose(np.atleast_2d(B_eq)), np.transpose(ub)), axis=0)

prob = osqp.OSQP()
prob.setup(P=P, A=A_eq_osqp, l=lb_osqp, u=ub_osqp, verbose=True, adaptive_rho_interval=50) #, adaptive_rho_fraction = 0)
res = prob.solve()
x_opt_osqp_prev = res.x
for i in range(100):
    prob = osqp.OSQP()
    prob.setup(P=P, A=A_eq_osqp, l=lb_osqp, u=ub_osqp, verbose=False, adaptive_rho_interval=50)
    res = prob.solve()
    x_opt_osqp = res.x
    if not(np.array_equal(x_opt_osqp, x_opt_osqp_prev)):
        aux = 1
    x_opt_osqp_prev = x_opt_osqp
    print(i)

aux= 1
