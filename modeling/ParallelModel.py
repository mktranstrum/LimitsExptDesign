import numpy as np
from BaseModel import BaseModel

class ParallelModel(BaseModel):

    def __init__(self, models):
        self.models = models
        self.Ms = [0] + [model.M for model in models]
        self.cumMs = np.cumsum(self.Ms)
        self.Nmodels = len(models)

        BaseModel.__init__(self, self.cumMs[-1], self.models[0].N)

        self.description = "Parallel Model of ("
        for model in self.models:
            self.description += " %s," %model.description
        self.description += ")"

    def r(self, x):

        ans = np.empty( ( self.M, ))

        for i, model in enumerate(self.models):
            ans[ self.cumMs[i]: self.cumMs[i+1]] = model.r(x)

        return ans

    def v(self, x, v):
        ans = np.empty( ( self.M, ) )
        for i, model in enumerate(self.models):
            ans[ self.cumMs[i]: self.cumMs[i+1]] = model.v(x, v)
        return ans
        
    def j(self, x):
        ans = np.empty( (self.M, self.N))
        for i, model in enumerate(self.models):
            ans[ self.cumMs[i]: self.cumMs[i+1],:] = model.j(x)
        return ans                                                                                        

    def Avv(self, x, v):
        ans = np.empty( ( self.M, ) )
        for i, model in enumerate(self.models):
            ans[ self.cumMs[i]: self.cumMs[i+1]] = model.Avv(x, v)
        return ans

    def Auv(self, x, u, v):
        ans = np.empty( ( self.M, ) )
        for i, model in enumerate(self.models):
            ans[ self.cumMs[i]: self.cumMs[i+1]] = model.Auv(x, u, v)
        return ans
        
        
