import numpy as np

# Finite Difference Estimtes

def v_FD(func, x, v, h, n = 1):
    if n == 1:
        return ( func(x + 0.5*h*v)  - func(x - 0.5*h*v) ) / h

def Avv_FD(func, x, v, h, n = 1):
    if n == 1:
        return ( func(x + h*v) + func(x - h*v) - 2.0*func(x) )/(h*h)

def Auv_FD(func, x, u, v, h, n = 1):
    if n == 1:
        return ( func(x + h*u + h*v) - func(x + h*u - h*v) - func(x - h*u + h*v) + func(x - h*u - h*v) )/(4*h*h)
        
class BaseModel(object):

    def __init__(self, M, N, description = "Base Model"):
        self.N = N
        self.M = M
        self.description = description
        self.h1 = np.sqrt( np.finfo(np.float64).eps ) # step size to estimate derivative
        self.h2 = np.sqrt( self.h1 ) # step size to estimate second derivatives
        self.n_v = 1 ## order of FD derivative estimate
        self.n_vv = 1 ## order of FD second directional second derivative estimate
        self.n_uv = 1 ## order of FD second derivative estimate
        
    def r(self, x):
        return x

    def v(self, x, v):
        return v_FD(self.r, x, v, self.h1, self.n_v)

    def Checkv_FD(self, x, v):
        return self.v(x,v) - v_FD(self.r, x, v, self.h1, self.n_v)

    def j(self, x):
        return np.array([ self.v(x, v) for v in np.eye(self.N) ] ).T

    def Avv(self, x, v):
        return Avv_FD(self.r, x, v, self.h2, self.n_vv)

    def CheckAvv_FD(self, x, v):
        return self.Avv(x,v) - Avv_FD(self.r, x, v, self.h2, self.n_vv)

    def Auv(self, x, u, v):
        return Auv_FD(self.r, x, u, v, self.h2, self.n_uv)

    def CheckAuv_FD(self, x, u, v):
        return self.Auv(x,u,v) - Auv_FD(self.r, x, u, v, self.h2, self.n_uv)

    def A(self, x):
        return np.array([[self.Auv(x, u, v) for u in np.eye(self.N)] for v in np.eye(self.N)]).T

        
