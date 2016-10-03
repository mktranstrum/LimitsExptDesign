import sys
import numpy as np

###########################
#Seeding the Random Numbers
###########################
seed = int( sys.argv[1] )
np.random.seed(seed)

import numpy as np

from modeling import BaseModel, DAECalc, Composition, ParallelModel

calc = DAECalc.DAECalc("__PC12_48__","__dPC12_48__","__d2PC12_48__")


class PC12_48(BaseModel):
    def __init__(self, ts, data = 0, weights = 1.0):
        self.ts = ts
        self.weights = weights
        self.data = data
        self.calc = calc #DAECalc.DAECalc("__PC12_48__","__dPC12_48__","__d2PC12_48__")
        BaseModel.__init__(self,len(ts) * 15, 69, "PC12_48") #15 dVars, 69 Parameters
        self.calc.kwargs['max_steps'] = 25000
        
    def r(self, x): #residual function
        return (self.calc.evaluate(x, self.ts).flatten() - self.data)*self.weights #A matrix flattened into a 300 value vector
        
    def v(self, x, v):
        return self.calc.evaluate_derivative(x, v, self.ts).flatten()*self.weights
        


class Expt(BaseModel):
    def __init__(self, x0):
        self.x0 = x0 # x0 contains the default values for the experiment
        BaseModel.__init__(self, 69, 48, "Expt")
        
    def r(self, x): #residuals
        return np.append(self.x0, x)
        
class Prior(BaseModel):
    def __init__(self, x0 = 1, weights = 25.0):
	self.x0 = x0
	self.weights = weights 
	BaseModel.__init__(self,48,48, "LinearPrior")

    def r(self, x):
	return np.log(x/(self.x0))*self.weights #(x - self.x0) * self.weights

    def v(self, x, v):
	return (v/x) * self.weights


#########################
### Setup base model
##########################
ts = np.linspace(0, 120, 100)
basemodel = PC12_48(ts)

#########################
###  First Apgar Expt
######################### 
expt1 = np.ones(21)
expt1[-1] = 0 # Zero
expt1[5] = 125.025*1e-2 # EGF
expt1[6] = 45.6*1e2 # NGF
expt1[9] *= 100 # Sos
expt1[11] *= 100 # Ras
expt1[18] *= 100 # C3G
Expt1Model = Expt(expt1)
basemodel1 = PC12_48(ts)
model1 = Composition(basemodel1, Expt1Model)

#########################
###  Second Apgar Expt
######################### 
expt2 = np.ones(21)
expt2[-1] = 0 # Zero
expt2[5] = 125.025*1e-6 # EGF
expt2[6] = 45.6*1e-4 # NGF
expt2[14] *= 100 # Mek
expt2[15] *= 100 # Erk
expt2[4] *= 0 # Raf1PPtase
Expt2Model = Expt(expt2)
basemodel2 = PC12_48(ts)
model2 = Composition(basemodel2, Expt2Model)


#########################
###  Third Apgar Expt
#########################
expt3 = np.ones(21)
expt3[-1] = 0 # Zero
expt3[5] = 0.0 # EGF
expt3[6] = 45.6*1e0 # NGF
expt3[13] *= 100 # BRaf
expt3[19] *= 100 # Rap1
expt3[2] *= 0 # RapGap
Expt3Model = Expt(expt3)
basemodel3 = PC12_48(ts)
model3 = Composition(basemodel3, Expt3Model)

#########################
###  Fourth Apgar Expt
#########################
expt4 = np.ones(21)
expt4[-1] = 0 # Zero
expt4[5] = 125.025*1e-6 # EGF
expt4[6] = 45.6*1e2 # NGF
expt4[10] *= 100 # P90Rsk
expt4[16] *= 100 # PI3K
expt4[17] *= 100 # Akt
Expt4Model = Expt(expt4)
basemodel4 = PC12_48(ts)
model4 = Composition(basemodel4, Expt4Model)

#########################
###  Fifth Apgar Expt
#########################
expt5 = np.ones(21)
expt5[-1] = 0 # Zero
expt5[5] = 125.025*1e-4 # EGF
expt5[6] = 45.6*1e-2 # NGF
expt5[12] *= 100 # Raf1
expt5[1] *= 0 # RasGap
Expt5Model = Expt(expt5)
basemodel5 = PC12_48(ts)
model5 = Composition(basemodel5, Expt5Model)

x = np.ones(48) #np.log(np.loadtxt("TrueParameters.txt")) 
#expx = np.exp(x) 

r = model1.r(x) #expx)
y1 = r.reshape(len(ts), 15)
sigma1 = np.maximum(0.1*y1.flatten(), 1e-2)  ## 10 percent error bar
basemodel1.weights = 1.0/np.loadtxt("apgardata/sigma1t.txt").flatten()
basemodel1.data = np.loadtxt("apgardata/y1t.txt").flatten()

r = model2.r(x) #expx)
y2 = r.reshape(len(ts), 15)
sigma2 = np.maximum(0.1*y2.flatten(), 1e-2)
basemodel2.weights = 1.0/np.loadtxt("apgardata/sigma2t.txt").flatten()
basemodel2.data = np.loadtxt("apgardata/y2t.txt").flatten()
 
r = model3.r(x) #expx)
y3 = r.reshape(len(ts), 15)
sigma3 = np.maximum(0.1*y3.flatten(), 1e-2)
basemodel3.weights = 1.0/np.loadtxt("apgardata/sigma3t.txt").flatten()
basemodel3.data = np.loadtxt("apgardata/y3t.txt").flatten()

r = model4.r(x) #expx)
y4 = r.reshape(len(ts), 15)
sigma4 = np.maximum(0.1*y4.flatten(), 1e-2)
basemodel4.weights = 1.0/np.loadtxt("apgardata/sigma4t.txt").flatten()
basemodel4.data = np.loadtxt("apgardata/y4t.txt").flatten()
 
r = model5.r(x) #expx)
y5 = r.reshape(len(ts), 15)
sigma5 = np.maximum(0.1*y5.flatten(), 1e-2)
basemodel5.weights = 1.0/np.loadtxt("apgardata/sigma5t.txt").flatten()
basemodel5.data = np.loadtxt("apgardata/y5t.txt").flatten()

PriorModel = Prior(weights = 1.0)

from Exponential import Exponential
model = Composition(ParallelModel([model1, model2, model3, model4, model5]), Exponential(48))
#x = np.log(np.loadtxt("TrueParameters.txt"))


from geodesiclm import geodesiclm

# Choose starting point
x0 = np.random.randn(48) * 1.0
PriorModel.x0 = np.exp(x0.copy())

xf, info = geodesiclm(model.r, x0, jacobian = model.j, full_output = 1, ibold = 0, print_level = 2, ibroyden = 1, artol = 0.01)

np.savetxt("fits/xf48_%03i.txt" %seed, xf) # Save the best fit
