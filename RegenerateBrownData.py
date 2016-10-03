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
        
#expts = np.loadtxt("Expts.txt")

# # expt = expts[0]
# expt = np.ones(21)
# expt[5] = 0.0 #125.025 # EGF
# expt[6] = 45.6
# expt[-1] = 0 # Zero

#########################
### Setup base model
##########################
ts = np.array([0,2,5,10,15,30,60,90,120])
basemodel = PC12_48(ts)

#########################
###  First Brown Expt
######################### 
# Wild type - only EGF stimulation
expt1 = np.ones(21)
expt1[-1] = 0 # Zero
expt1[5] = 125.025 # EGF 
expt1[6] = 0 # NGF 
BrownModel1 = Expt(expt1)
basemodel1 = PC12_48(ts)
model1 = Composition(basemodel1, BrownModel1)

#########################
###  Second Brown Expt
######################### 
# Wild type - only NGF stimulation
expt2 = np.ones(21)
expt2[-1] = 0 # Zero
expt2[5] = 0 # EGF 
expt2[6] = 45.6 # NGF 
BrownModel2 = Expt(expt2)
basemodel2 = PC12_48(ts)
model2 = Composition(basemodel2, BrownModel2)


#########################
###  Third Brown Expt
#########################
# EGF stimulation, PI3K knockout
expt3 = np.ones(21)
expt3[-1] = 0 # Zero
expt3[5] = 125.025 # EGF 
expt3[6] = 0 # NGF 
expt3[16] *= 0 # Total_PI3K (links to PI3K from Ras and EGFR cut)
BrownModel3 = Expt(expt3)
basemodel3 = PC12_48(ts)
model3 = Composition(basemodel3, BrownModel3)


#########################
###  Fourth Brown Expt
######################### 
# NGF stimulation, PI3K knockout
expt4 = np.ones(21)
expt4[-1] = 0 # Zero
expt4[5] = 0 # EGF 
expt4[6] = 45.6 # NGF 
expt4[16] *= 0 # Total_PI3K (links to PI3K from Ras)
BrownModel4 = Expt(expt4)
basemodel4 = PC12_48(ts)
model4 = Composition(basemodel4, BrownModel4)

#########################
###  Fifth Brown Expt
#########################
# EGF stimulation, Raf1 knockout
#Won't this just give zero Erk?
expt5 = np.ones(21)
expt5[-1] = 0 # Zero
expt5[5] = 125.025 # EGF 
expt5[6] = 0 # NGF 
expt5[12] *= 0 # Total_Raf1 (links to PI3K from Ras)
BrownModel5 = Expt(expt5)
basemodel5 = PC12_48(ts)
model5 = Composition(basemodel5, BrownModel5)

#########################
###  Sixth Brown Expt
######################### 
# NGF stimulation, Raf1 knockout
expt6 = np.ones(21)
expt6[-1] = 0 # Zero
expt6[5] = 0 # EGF 
expt6[6] = 45.6 # NGF 
expt6[12] *= 0 # Total_Raf1 (links to Raf1 from Ras)
BrownModel6 = Expt(expt6)
basemodel6 = PC12_48(ts)
model6 = Composition(basemodel6, BrownModel6)

#########################
###  Seventh Brown Expt
#########################
# NGF stimulation, B-Raf knockout
expt7 = np.ones(21)
expt7[-1] = 0 # Zero
expt7[5] = 0 # EGF 
expt7[6] = 45.6 # NGF 
expt7[13] *= 0 # Total_BRaf (links to B-Raf from Rap1)
BrownModel7 = Expt(expt7)
basemodel7 = PC12_48(ts)
model7 = Composition(basemodel7, BrownModel7)

x = np.log(np.loadtxt("TrueParameters.txt")) 
expx = np.exp(x) 

r = model1.r(expx)
y1 = r.reshape(len(ts), 15)
sigma1 = np.maximum(0.1*y1, 1e-2)  ## 10 percent error bar
weight1 = 1.0/sigma1
weight1[:,0] *= 0 # weights for EGF -> 0
weight1[:,1] *= 0 #weights for NGF -> 0
weight1[:,2] *= 0 #weights for boundEGFReceptor -> 0
weight1[:,3] *= 0 #weights for boundNGFReceptor -> 0
weight1[:,4] *= 0 #weights for SosActive -> 0
weight1[:,5] *= 0 #weights for P90RskActive -> 0
weight1[:,11] *= 0 #weights for PI3KActive -> 0
weight1[:,12] *= 0 #weights for AktActive -> 0
weight1[:,13] *= 0 #weights for CSGActive -> 0
# repeat for all other non-observed variables
basemodel1.weights = weight1.copy().flatten()
np.savetxt("browndata/y1.txt", y1[:,:15]) # Noise: +sigma1... 
np.savetxt("browndata/weight1.txt", weight1[:,:15])


r = model2.r(expx)
y2 = r.reshape(len(ts), 15)
sigma2 = np.maximum(0.1*y1, 1e-2)  ## 10 percent error bar
weight2 = 1.0/sigma2
weight2[:,0] *= 0 # weights for EGF -> 0
weight2[:,1] *= 0 #weights for NGF -> 0
weight2[:,2] *= 0 #weights for boundEGFReceptor -> 0
weight2[:,3] *= 0 #weights for boundNGFReceptor -> 0
weight2[:,4] *= 0 #weights for SosActive -> 0
weight2[:,5] *= 0 #weights for P90RskActive -> 0
weight2[:,11] *= 0 #weights for PI3KActive -> 0
weight2[:,12] *= 0 #weights for AktActive -> 0
weight2[:,13] *= 0 #weights for CSGActive -> 0
# repeat for all other non-observed variables
basemodel2.weights = weight2.copy().flatten()
np.savetxt("browndata/y2.txt", y2[:,:15]) # Noise: +sigma1... 
np.savetxt("browndata/weight2.txt", weight2[:,:15])

 
r = model3.r(expx)
y3 = r.reshape(len(ts), 15)
sigma3 = np.maximum(0.1*y3, 1e-2)
weight3 = 1.0/sigma3
weight3[:,0] *= 0 # weights for EGF -> 0
weight3[:,1] *= 0 #weights for NGF -> 0
weight3[:,2] *= 0 #weights for boundEGFReceptor -> 0
weight3[:,3] *= 0 #weights for boundNGFReceptor -> 0
weight3[:,4] *= 0 #weights for SosActive -> 0
weight3[:,5] *= 0 #weights for P90RskActive -> 0
weight3[:,11] *= 0 #weights for PI3KActive -> 0
weight3[:,12] *= 0 #weights for AktActive -> 0
weight3[:,13] *= 0 #weights for CSGActive -> 0
# repeat for all other non-observed variables
basemodel3.weights = weight3.copy().flatten()
np.savetxt("browndata/y3.txt", y3[:,:15]) # Noise: +sigma1... 
np.savetxt("browndata/weight3.txt", weight3[:,:15])


r = model4.r(expx)
y4 = r.reshape(len(ts), 15)
sigma4 = np.maximum(0.1*y4, 1e-2)
weight4 = 1.0/sigma4
weight4[:,0] *= 0 # weights for EGF -> 0
weight4[:,1] *= 0 #weights for NGF -> 0
weight4[:,2] *= 0 #weights for boundEGFReceptor -> 0
weight4[:,3] *= 0 #weights for boundNGFReceptor -> 0
weight4[:,4] *= 0 #weights for SosActive -> 0
weight4[:,5] *= 0 #weights for P90RskActive -> 0
weight4[:,11] *= 0 #weights for PI3KActive -> 0
weight4[:,12] *= 0 #weights for AktActive -> 0
weight4[:,13] *= 0 #weights for CSGActive -> 0
# repeat for all other non-observed variables
basemodel4.weights = weight4.copy().flatten()
np.savetxt("browndata/y4.txt", y4[:,:15]) # Noise: +sigma1... 
np.savetxt("browndata/weight4.txt", weight4[:,:15])

 
r = model5.r(expx)
y5 = r.reshape(len(ts), 15)
sigma5 = np.maximum(0.1*y5, 1e-2)
weight5 = 1.0/sigma5
weight5[:,0] *= 0 # weights for EGF -> 0
weight5[:,1] *= 0 #weights for NGF -> 0
weight5[:,2] *= 0 #weights for boundEGFReceptor -> 0
weight5[:,3] *= 0 #weights for boundNGFReceptor -> 0
weight5[:,4] *= 0 #weights for SosActive -> 0
weight5[:,5] *= 0 #weights for P90RskActive -> 0
weight5[:,11] *= 0 #weights for PI3KActive -> 0
weight5[:,12] *= 0 #weights for AktActive -> 0
weight5[:,13] *= 0 #weights for CSGActive -> 0
# repeat for all other non-observed variables
basemodel5.weights = weight5.copy().flatten()
np.savetxt("browndata/y5.txt", y5[:,:15]) # Noise: +sigma1... 
np.savetxt("browndata/weight5.txt", weight5[:,:15])


r = model6.r(expx)
y6 = r.reshape(len(ts), 15)
sigma6 = np.maximum(0.1*y1, 1e-2)  ## 10 percent error bar
weight6 = 1.0/sigma6
weight6[:,0] *= 0 # weights for EGF -> 0
weight6[:,1] *= 0 #weights for NGF -> 0
weight6[:,2] *= 0 #weights for boundEGFReceptor -> 0
weight6[:,3] *= 0 #weights for boundNGFReceptor -> 0
weight6[:,4] *= 0 #weights for SosActive -> 0
weight6[:,5] *= 0 #weights for P90RskActive -> 0
weight6[:,11] *= 0 #weights for PI3KActive -> 0
weight6[:,12] *= 0 #weights for AktActive -> 0
weight6[:,13] *= 0 #weights for CSGActive -> 0
# repeat for all other non-observed variables
basemodel6.weights = weight6.copy().flatten()
#np.savetxt("browndata/y1.txt", y1[:15]+sigma1... 
np.savetxt("browndata/y6.txt", y6[:,:15]) # Noise: +sigma1... 
np.savetxt("browndata/weight6.txt", weight6[:,:15])


r = model7.r(expx)
y7 = r.reshape(len(ts), 15)
sigma7 = np.maximum(0.1*y1, 1e-2)  ## 10 percent error bar
weight7 = 1.0/sigma7
weight7[:,0] *= 0 # weights for EGF -> 0
weight7[:,1] *= 0 #weights for NGF -> 0
weight7[:,2] *= 0 #weights for boundEGFReceptor -> 0
weight7[:,3] *= 0 #weights for boundNGFReceptor -> 0
weight7[:,4] *= 0 #weights for SosActive -> 0
weight7[:,5] *= 0 #weights for P90RskActive -> 0
weight7[:,11] *= 0 #weights for PI3KActive -> 0
weight7[:,12] *= 0 #weights for AktActive -> 0
weight7[:,13] *= 0 #weights for CSGActive -> 0
# repeat for all other non-observed variables
basemodel7.weights = weight7.copy().flatten()
np.savetxt("browndata/y7.txt", y7[:,:15]) # Noise: +sigma1... 
np.savetxt("browndata/weight7.txt", weight7[:,:15])


from Exponential import Exponential
model = Composition(ParallelModel([model1, model2, model3, model4, model5, model6, model7]), Exponential(48))


