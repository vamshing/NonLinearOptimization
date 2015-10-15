
# coding: utf-8

# In[9]:

import math
import numpy as np
from sympy import *


# In[10]:

x0 = np.array([1,1,1,1])


# In[11]:

#function Defition

def func(a,b,c,d):
    x,y,z,w=symbols('x y z w')
    f = (x - 2*y**2)**2 + (y - 2*z**2)**2 + (z - 2*w**2)**2
    return f.subs({x:a,y:b,z:c,w:d})

def derv_x(a,b,c,d):
    x,y,z,w=symbols('x y z w')
    f = (x - 2*y**2)**2 + (y - 2*z**2)**2 + (z - 2*w**2)**2
    return diff(f,x).subs({x:a,y:b,z:c,w:d})
    
def derv_y(a,b,c,d):
    x,y,z,w=symbols('x y z w')
    f = (x - 2*y**2)**2 + (y - 2*z**2)**2 + (z - 2*w**2)**2
    return diff(f,y).subs({x:a,y:b,z:c,w:d})    
    
def derv_z(a,b,c,d):
    x,y,z,w=symbols('x y z w')
    f = (x - 2*y**2)**2 + (y - 2*z**2)**2 + (z - 2*w**2)**2
    return diff(f,z).subs({x:a,y:b,z:c,w:d})
    
def derv_w(a,b,c,d):
    x,y,z,w=symbols('x y z w')
    f = (x - 2*y**2)**2 + (y - 2*z**2)**2 + (z - 2*w**2)**2
    return diff(f,w).subs({x:a,y:b,z:c,w:d})

def grad_vector(x0):
    return np.matrix([derv_x(x0[0],x0[1],x0[2],x0[3]),derv_y(x0[0],x0[1],x0[2],x0[3]),
                         derv_z(x0[0],x0[1],x0[2],x0[3]),derv_w(x0[0],x0[1],x0[2],x0[3])])


# In[12]:

grad_vector(x0)


# In[42]:

def non_linear_conjugate_gradient(A, b, x0, tol = 1.0e-8, max_iter = 4):
    """
    A function to solve [A]{x} = {b} linear equation system with the 
    conjugate gradient method.
    
    :param A : array 
        A real symmetric positive definite matrix(assumed)
        
    :param b : vector
        The vector of the system which is given in RHS.
        
    :param x0 : vector
        The starting guess for the solution.
        
    :param max_iter : integer
        Maximum number of iterations. Iteration will stop after max_iter 
        steps even if the specified tolerance has not been achieved.
        
    :param tol : float
        Tolerance to achieve. The algorithm will terminate when either 
        the relative or the absolute residual is below tol.
        
    :var    r0 : vector
                 Initialization stores the value (b - a * A )
    
    :var    d  : vector
    
    :var    a  : float
                 Iteratively computes the scalar of (r1T.r1)/(r0T.r0)

    :var    ri : vector
                 Iteratively stores the value (r - a * A * d), used to check for the convergence
    
    :var    x  : vector 
                 Stores the solution for the next iteration iteratively
                 
    """
    x = x0
    r0 = grad_vector(x0) * (-1)
    d = r0
    alpha = 0.5

#   Iterations:   
    for i in xrange(max_iter):
        
        #compute the lamda value by line search
        lmda = Symbol('lmda')
    # function value at x0
        f_x0 = func(x0)
    # The Goldstien-Armijo criteria for the lambda selection
        rhs = f_x0 + np.dot((d.T)*(alpha)*(lmda),grad_vector(x0))
        lhs = f.subs({x:guess.item(0)+lmda*d.item(0),y:guess.item(1)+lmda*d.item(1)}) 
    # solver for the lamda value from quadratic inequality
        try:
            lmda_value = max(solve(lhs-rhs,lmda))
        except ValueError: 
            pass
        
        # line search method ends there
        x = x + d*lmda
        
        print "iteration: ",i, "r(i): ",round(np.linalg.norm(grad_vector(x)),5)

        if np.linalg.norm(grad_vector(x)) < tol:
            print "\nConverged Successfully in iterations :",i
            print "The result of vector x:"
            return np.around(x,decimals=10)
            break
        b = float((grad_vector(x)*grad_vector(x).T)/(grad_vector(x0)*grad_vector(x0).T))
        d = grad_vector(x) * (-1) + b * d
        x0 = x
    return x


# In[62]:

x0 = [1,1,1,1]
x = x0
r0 = grad_vector(x0) * (-1)
d = r0
alpha = 0.5


# In[64]:

lmda = Symbol('lmda')
f_x0 = func(x0[0],x0[1],x0[2],x0[3])
rhs = f_x0 + 


# In[66]:

np.dot((d.T)*(alpha)*(lmda),grad_vector(x0))


# In[74]:

d.T * alpha * lmda


# In[ ]:



