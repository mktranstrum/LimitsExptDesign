
import numpy as np
import sys
from modeling import BaseModel, DAECalc, Composition, ParallelModel

###########################
#Seeding the Random Numbers
###########################
seed = int( sys.argv[1] )
np.random.seed(seed)



calc = DAECalc.DAECalc("__PC12_MA__","__dPC12_MA__","__d2PC12_MA__")


class PC12_MA(BaseModel):
    def __init__(self, ts, data = 0,  weights = 1.0):
        self.ts = ts
        self.weights = weights
	self.data = data
        self.calc = calc #DAECalc.DAECalc("__PC12_MA__","__dPC12_MA__","__d2PC12_MA__")
        BaseModel.__init__(self,len(ts) * 15, 91, "PC12_MA") #15 of 54 dVars, 91 Parameters
        self.calc.kwargs['max_steps']=5000
        
    def r(self, x): #residual function
        return (((self.calc.evaluate(x, self.ts)[:,:15]-self.data)*self.weights).flatten()) #A matrix flattened into a 1500 value vector
        
    def v(self, x, v):
        return (self.calc.evaluate_derivative(x, v, self.ts)[:,:15]*self.weights).flatten()
        


class Expt(BaseModel):
    def __init__(self, x0):
        self.x0 = x0 # x0 contains the default values for the experiment
        BaseModel.__init__(self, 91, 70, "Expt")
        
    def r(self, x): #residuals
        return np.append(self.x0, x)
        
class Prior(BaseModel):
    def __init__(self, x0 = 1, weights = 25.0):
	self.x0 = x0
	self.weights = weights 
	BaseModel.__init__(self,70,70, "LinearPrior")

    def r(self, x):
	return np.log(x/(self.x0))*self.weights #(x - self.x0) * self.weights

    def v(self, x, v):
	return (v/x) * self.weights

class MMPrior(BaseModel):
    def __init__(self, weights):
	self.weights  = weights
	BaseModel.__init__(self, 22, 70, "MMPrior")

    def r(self,x):
	alpha = 10.0
#    General Form:
#	r1 = (kcat/kr) * self.weights #Equilibrium Approximation
#	r2 = (kf / (kr + kcat)) * self.weights #QSSA
#	r = alpha*(tanh(r1/alpha)) * alpha*(tanh(r2/alpha))
#    Specific cases:
		
	d = 21
	k02f = x[25-d]
	k02b = x[26-d]
	k03f = x[27-d]

	k04f = x[28-d]
	k04b = x[29-d]
	k05f = x[30-d]

	k06f = x[31-d]
	k06b = x[32-d]
	k07f = x[33-d]

	k08f = x[34-d]
	k08b = x[35-d]
	k09f = x[36-d]
	
	k10f = x[37-d]
	k10b = x[38-d]
	k11f = x[39-d]

	k12f = x[40-d]
	k12b = x[41-d]
	k13f = x[42-d]

	k14f = x[43-d]
	k14b = x[44-d]
	k15f = x[45-d]

	k16f = x[46-d]
	k16b = x[47-d]
	k17f = x[48-d]

	k18f = x[49-d]
	k18b = x[50-d]
	k19f = x[51-d]

	k20f = x[52-d]
	k20b = x[53-d]
	k21f = x[54-d]

	k22f = x[55-d]
	k22b = x[56-d]
	k23f = x[57-d]

	k24f = x[58-d]
	k24b = x[59-d]
	k25f = x[60-d]

	k26f = x[61-d]
	k26b = x[62-d]
	k27f = x[63-d]

	k28f = x[64-d]
	k28b = x[65-d]
	k29f = x[66-d]

	k30f = x[67-d]
	k30b = x[68-d]
	k31f = x[69-d]

	k32f = x[70-d]
	k32b = x[71-d]
	k33f = x[72-d]

	k34f = x[73-d]
	k34b = x[74-d]
	k35f = x[75-d]

	k36f = x[76-d]
	k36b = x[77-d]
	k37f = x[78-d]

	k38f = x[79-d]
	k38b = x[80-d]
	k39f = x[81-d]

	k40f = x[82-d]
	k40b = x[83-d]
	k41f = x[84-d]

	k42f = x[85-d]
	k42b = x[86-d]
	k43f = x[87-d]

	k44f = x[88-d]
	k44b = x[89-d]
	k45f = x[90-d]

	rs = []

	rE1 = (k03f/k02b) * self.weights
	rQ1 = (k02f/(k02b + k03f)) * self.weights 
	r01 = alpha*np.tanh(rE1/alpha) * alpha*np.tanh(rQ1/alpha)
	rs.append(r01)	

	rE2 = (k05f/k04b) * self.weights
	rQ2 = (k04f/(k04b + k05f)) * self.weights 
	r02 = alpha*np.tanh(rE2/alpha) * alpha*np.tanh(rQ2/alpha)
	rs.append(r02)	

	rE3 = (k07f/k06b) * self.weights
	rQ3 = (k06f/(k06b + k07f)) * self.weights 
	r03 = alpha*np.tanh(rE3/alpha) * alpha*np.tanh(rQ3/alpha)
	rs.append(r03)	

	rE4 = (k09f/k08b) * self.weights
	rQ4 = (k08f/(k08b + k09f)) * self.weights 
	r04 = alpha*np.tanh(rE4/alpha) * alpha*np.tanh(rQ4/alpha)
	rs.append(r04)	

	rE5 = (k11f/k10b) * self.weights
	rQ5 = (k10f/(k10b + k11f)) * self.weights 
	r05 = alpha*np.tanh(rE5/alpha) * alpha*np.tanh(rQ5/alpha)
	rs.append(r05)	

	rE6 = (k13f/k12b) * self.weights
	rQ6 = (k12f/(k12b + k13f)) * self.weights 
	r06 = alpha*np.tanh(rE6/alpha) * alpha*np.tanh(rQ6/alpha)
	rs.append(r06)	

	rE7 = (k15f/k14b) * self.weights
	rQ7 = (k14f/(k14b + k15f)) * self.weights 
	r07 = alpha*np.tanh(rE7/alpha) * alpha*np.tanh(rQ7/alpha)
	rs.append(r07)	

	rE8 = (k17f/k16b) * self.weights
	rQ8 = (k16f/(k16b + k17f)) * self.weights 
	r08 = alpha*np.tanh(rE8/alpha) * alpha*np.tanh(rQ8/alpha)
	rs.append(r08)	

	rE9 = (k19f/k18b) * self.weights
	rQ9 = (k18f/(k18b + k19f)) * self.weights 
	r09 = alpha*np.tanh(rE9/alpha) * alpha*np.tanh(rQ9/alpha)
	rs.append(r09)	

	rE10 = (k21f/k20b) * self.weights
	rQ10 = (k20f/(k20b + k21f)) * self.weights 
	r10 = alpha*np.tanh(rE10/alpha) * alpha*np.tanh(rQ10/alpha)
	rs.append(r10)	

	rE11 = (k23f/k22b) * self.weights
	rQ11 = (k22f/(k22b + k23f)) * self.weights 
	r11 = alpha*np.tanh(rE11/alpha) * alpha*np.tanh(rQ11/alpha)
	rs.append(r11)	

	rE12 = (k25f/k24b) * self.weights
	rQ12 = (k24f/(k24b + k25f)) * self.weights 
	r12 = alpha*np.tanh(rE12/alpha) * alpha*np.tanh(rQ12/alpha)
	rs.append(r12)	

	rE13 = (k27f/k26b) * self.weights
	rQ13 = (k26f/(k26b + k27f)) * self.weights 
	r13 = alpha*np.tanh(rE13/alpha) * alpha*np.tanh(rQ13/alpha)
	rs.append(r13)	

	rE14 = (k29f/k28b) * self.weights
	rQ14 = (k28f/(k28b + k29f)) * self.weights 
	r14 = alpha*np.tanh(rE14/alpha) * alpha*np.tanh(rQ14/alpha)
	rs.append(r14)	

	rE15 = (k31f/k30b) * self.weights
	rQ15 = (k30f/(k30b + k31f)) * self.weights 
	r15 = alpha*np.tanh(rE15/alpha) * alpha*np.tanh(rQ15/alpha)
	rs.append(r15)	

	rE16 = (k33f/k32b) * self.weights
	rQ16 = (k32f/(k32b + k33f)) * self.weights 
	r16 = alpha*np.tanh(rE16/alpha) * alpha*np.tanh(rQ16/alpha)
	rs.append(r16)	

	rE17 = (k35f/k34b) * self.weights
	rQ17 = (k34f/(k34b + k35f)) * self.weights 
	r17 = alpha*np.tanh(rE17/alpha) * alpha*np.tanh(rQ17/alpha)
	rs.append(r17)	

	rE18 = (k37f/k36b) * self.weights
	rQ18 = (k36f/(k36b + k37f)) * self.weights 
	r18 = alpha*np.tanh(rE18/alpha) * alpha*np.tanh(rQ18/alpha)
	rs.append(r18)	

	rE19 = (k39f/k38b) * self.weights
	rQ19 = (k38f/(k38b + k39f)) * self.weights 
	r19 = alpha*np.tanh(rE19/alpha) * alpha*np.tanh(rQ19/alpha)
	rs.append(r19)	

	rE20 = (k41f/k40b) * self.weights
	rQ20 = (k40f/(k40b + k41f)) * self.weights 
	r20 = alpha*np.tanh(rE20/alpha) * alpha*np.tanh(rQ20/alpha)
	rs.append(r20)	

	rE21 = (k43f/k42b) * self.weights
	rQ21 = (k42f/(k42b + k43f)) * self.weights 
	r21 = alpha*np.tanh(rE21/alpha) * alpha*np.tanh(rQ21/alpha)
	rs.append(r21)	

	rE22 = (k45f/k44b) * self.weights
	rQ22 = (k44f/(k44b + k45f)) * self.weights 
	r22 = alpha*np.tanh(rE22/alpha) * alpha*np.tanh(rQ22/alpha)
	rs.append(r22)	

	return np.array(rs)

#########################
### Setup base model
##########################
ts = np.array([0,2,5,10,15,30,60,90,120]) #np.linspace(0, 120, 100)
basemodel = PC12_MA(ts)

#########################
###  First Brown Expt
######################### 
# Wild type - only EGF stimulation
expt1 = np.ones(21)
expt1[-1] = 0 # Zero
expt1[5] = 125.025 # EGF 
expt1[6] = 0 # NGF 
BrownModel1 = Expt(expt1)
sigma = 1.0/(np.loadtxt("browndata/weight1.txt")+ 1e-10)
basemodel1 = PC12_MA(ts, np.loadtxt("browndata/y1.txt")+ np.random.randn(9,15)*sigma, np.loadtxt("browndata/weight1.txt"))
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
sigma = 1.0/(np.loadtxt("browndata/weight2.txt")+ 1e-10)
basemodel2 = PC12_MA(ts, np.loadtxt("browndata/y2.txt")+np.random.randn(9,15)*sigma, np.loadtxt("browndata/weight2.txt"))
#basemodel2 = PC12_MA(ts)
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
sigma = 1.0/(np.loadtxt("browndata/weight3.txt")+ 1e-10)
basemodel3 = PC12_MA(ts, np.loadtxt("browndata/y3.txt") + np.random.randn(9,15)*sigma, np.loadtxt("browndata/weight3.txt"))
#basemodel3 = PC12_MA(ts)
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
sigma = 1.0/(np.loadtxt("browndata/weight4.txt")+ 1e-10)
basemodel4 = PC12_MA(ts, np.loadtxt("browndata/y4.txt")+np.random.randn(9,15)*sigma, np.loadtxt("browndata/weight4.txt"))
#basemodel4 = PC12_MA(ts)
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
sigma = 1.0/(np.loadtxt("browndata/weight5.txt")+ 1e-10)
basemodel5 = PC12_MA(ts, np.loadtxt("browndata/y5.txt")+np.random.randn(9,15)*sigma, np.loadtxt("browndata/weight5.txt"))
#basemodel5 = PC12_MA(ts)
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
sigma = 1.0/(np.loadtxt("browndata/weight6.txt")+ 1e-10)
basemodel6 = PC12_MA(ts, np.loadtxt("browndata/y6.txt")+np.random.randn(9,15)*sigma, np.loadtxt("browndata/weight6.txt"))
#basemodel6 = PC12_MA(ts)
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
sigma = 1.0/(np.loadtxt("browndata/weight7.txt")+ 1e-10)
basemodel7 = PC12_MA(ts, np.loadtxt("browndata/y7.txt")+np.random.randn(9,15)*sigma, np.loadtxt("browndata/weight7.txt"))
#basemodel7 = PC12_MA(ts)
model7 = Composition(basemodel7, BrownModel7)


MMPriorModel = MMPrior(10)  # Prior to force the Km parameter combinations to be large
# Regularizing prior to prevent running into boundaries prematuraly, inspired by: The geometry of nonlinear least squares with applications to sloppy models and optimization. Physical Review E. 83 036701, (2011)
PriorModel = Prior(1, 1) # Change these arguments to vary the location and strength of the regularizing term


#Currently: MMPrior included
from Exponential import Exponential
model = Composition(ParallelModel([model1, model2, model3, model4, model5, model6, model7, MMPriorModel, PriorModel]), Exponential(70)) 

from geodesiclm import geodesiclm

xi = np.random.randn(70) * 0.1 #  Choose a random starting point.
PriorModel.x0 = np.exp(xi)

xf, info = geodesiclm( model.r, xi, jacobian = model.j, full_output = 1, print_level = 2, ibold = 0, ibroyden = 1, artol = 1e-2)

np.savetxt("fits/xf_%03i_MMPrior.txt" %seed, xf) # Save the best fit


