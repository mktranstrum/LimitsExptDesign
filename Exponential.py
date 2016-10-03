import numpy as np
from modeling import BaseModel

class Exponential(BaseModel):
    """
    Exponential transformation of parameters
    """
    def __init__(self, N):
        BaseModel.__init__(self, N, N, "%i Parameter Exponential" %N)

    def r(self, x):
        return np.exp(x)

    def v(self, x, v):
        return np.exp(x)*v

    def Avv(self, x, v):
        return np.exp(x)*v*v

    def Auv(self, x, u, v):
        return np.exp(x)*u*v
