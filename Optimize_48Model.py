# -*- coding: utf-8 -*-
"""
Created on Thu Oct 09 11:05:42 2014

@author: Andrew
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Aug 19 19:28:35 2014

@author: Andrew
"""

import Model_PC12_48
reload(Model_PC12_48)
model = Model_PC12_48.model
import numpy as np
from geodesiclm import geodesiclm
import scipy

import sys
N = int(sys.argv[1])
seed = int(sys.argv[2])
np.random.seed(seed)
PriorModel = Model_PC12_48.PriorModel
x0 = np.random.randn(48) * 1.0
PriorModel.x0 = np.exp(x0.copy())

xf, info = geodesiclm(model.r, x0, jacobian = model.j, full_output = 1, ibold = 0, print_level = 2, ibroyden = 1, artol = 0.01)

np.savetxt("Fits/xf_%03i_%03i_%f_noMMPrior.txt" %(N, seed, PriorModel.weights), xf)
np.savetxt("Costs/Cost_%03i_%03i_%f_noMMPrior.txt" %(N, seed, PriorModel.weights), [sum(info['fvec']**2)/2] )
np.savetxt("Residuals/res_%03i_%03i_%f_noMMPrior.txt" %(N, seed, PriorModel.weights), info['fvec'] )

#def Optimize_Model(Model, n_parameters):
    
#    x0 = np.zeros(n_parameters)
    
#    x_bestfit = scipy.optimize.leastsq(Model.r, x0, Dfun=Model.j)

#    np.savetxt("48_BestFitNew.txt", x_bestfit)

#    return x_bestfit    

#best_fit = Optimize_Model(model, 48)

#print best_fit
