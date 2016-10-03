import numpy as np

# Edit this line to choose which parameter values to load
x = np.exp( np.loadtxt("fits/xf_%03i_MMPrior.txt" %0 )) 

# Edit this line to seed the random number generator for adding random noise to data
np.random.seed(0)


from modeling import BaseModel, DAECalc, Composition, ParallelModel

calc = DAECalc.DAECalc("__PC12_MA__","__dPC12_MA__","__d2PC12_MA__")


class PC12_MA(BaseModel):
    def __init__(self, ts, weights = 1.0):
        self.ts = ts
        self.weights = weights
        self.calc = calc #DAECalc.DAECalc("__PC12_MA__","__dPC12_MA__","__d2PC12_MA__")
        BaseModel.__init__(self,len(ts) * 54, 91, "PC12_MA") #54 (not 51) dVars, 91(21+70, not 21+64=85) Parameters
        self.calc.kwargs['max_steps']=5000
        
    def r(self, x): #residual function
        return self.calc.evaluate(x, self.ts).flatten()*self.weights #A matrix flattened into a 300 value vector
        
    def v(self, x, v):
        return self.calc.evaluate_derivative(x, v, self.ts).flatten()*self.weights
        


class Expt(BaseModel):
    def __init__(self, x0):
        self.x0 = x0 # x0 contains the default values for the experiment
        BaseModel.__init__(self, 91, 70, "Expt")
        
    def r(self, x): #residuals
        return np.append(self.x0, x)
        
#########################
### Setup base model
##########################
ts = np.linspace(0, 120, 100)
basemodel = PC12_MA(ts)

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
basemodel1 = PC12_MA(ts)
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
basemodel2 = PC12_MA(ts)
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
basemodel3 = PC12_MA(ts)
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
basemodel4 = PC12_MA(ts)
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
basemodel5 = PC12_MA(ts)
model5 = Composition(basemodel5, Expt5Model)

from Exponential import Exponential
model = ParallelModel([model1, model2, model3, model4, model5]) #Composed with Exponential 70 on 01/13/15

r = model1.r(x)
y1 = r.reshape(len(ts), 54)
sigma1 = np.maximum(0.1*y1, 1e-2)  ## 10 percent error bar
np.savetxt("apgardata/y1t.txt", y1[:,:15] + sigma1[:,:15]*np.random.randn(len(ts), 15) )
np.savetxt("apgardata/sigma1t.txt", sigma1[:,:15])
basemodel1.weights = 1.0/sigma1.flatten()

r = model2.r(x)
y2 = r.reshape(len(ts), 54)
sigma2 = np.maximum(0.1*y2, 1e-2)
np.savetxt("apgardata/y2t.txt", y2[:,:15] + sigma2[:,:15]*np.random.randn(len(ts),15) )
np.savetxt("apgardata/sigma2t.txt", sigma2[:,:15])
basemodel2.weights = 1.0/sigma2.flatten()

r = model3.r(x)
y3 = r.reshape(len(ts), 54)
sigma3 = np.maximum(0.1*y3, 1e-2)
np.savetxt("apgardata/y3t.txt", y3[:,:15] + sigma3[:,:15]*np.random.randn(len(ts),15) )
np.savetxt("apgardata/sigma3t.txt", sigma3[:,:15])
basemodel3.weights = 1.0/sigma3.flatten()

r = model4.r(x)
y4 = r.reshape(len(ts), 54)
sigma4 = np.maximum(0.1*y4, 1e-2)
np.savetxt("apgardata/y4t.txt", y4[:,:15] + sigma4[:,:15]*np.random.randn(len(ts),15) )
np.savetxt("apgardata/sigma4t.txt", sigma4[:,:15])
basemodel4.weights = 1.0/sigma4.flatten()

r = model5.r(x)
y5 = r.reshape(len(ts), 54)
sigma5 = np.maximum(0.1*y5, 1e-2)
np.savetxt("apgardata/y5t.txt", y5[:,:15] + sigma5[:,:15]*np.random.randn(len(ts),15) )
np.savetxt("apgardata/sigma5t.txt", sigma5[:,:15])
basemodel5.weights = 1.0/sigma5.flatten()


