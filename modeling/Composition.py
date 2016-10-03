import numpy as np
from BaseModel import BaseModel

class Composition(BaseModel):

    def __init__(self, modelout, modelin):
        BaseModel.__init__(self, modelout.M, modelin.N, "Composition of (%s) and (%s)" %(modelout.description, modelin.description) )
        self.modelout = modelout
        self.modelin = modelin
        ## How to calculate the jacobian
        ## calculatejac = True, means j = dot(jout, jin)
        ## calcualtejac = FAlse means j = array(v)
        ## If modelin.M > modelin.N, then the second is faster
        if modelin.M <= modelin.N:
            self.calculatejac = True
        else:
            self.calculatejac = False

    def r(self, x):
        return self.modelout.r(self.modelin.r(x))

    def v(self, x, v):
        return self.modelout.v(self.modelin.r(x), self.modelin.v(x, v) )

    def j(self, x):
        if self.calculatejac:
            return np.dot( self.modelout.j(self.modelin.r(x)), self.modelin.j(x) )
        else:
            return np.array([ self.v(x, v) for v in np.eye(self.N) ] ).T


    def Auv(self, x, u, v):
        return self.modelout.Auv( self.modelin.r(x), self.modelin.v(x, u), self.modelin.v(x,v) ) + self.modelout.v(self.modelin.r(x), self.modelin.Auv(x, u, v) )

    def Avv(self, x, v):
        return self.modelout.Avv(self.modelin.r(x), self.modelin.v(x,v) ) + self.modelout.v(self.modelin.r(x), self.modelin.Avv(x, v))
        
def MultipleComposition(*models):
    ans = Composition(models[-2], models[-1])
    for i in range(len(models) - 3, -1, -1):
        ans = Composition(models[i], ans)
    return ans
