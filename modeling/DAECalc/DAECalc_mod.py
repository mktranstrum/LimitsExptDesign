import sys, os, platform, importlib
import scipy
import scipy.linalg
from daskr import daeint
from scipy.optimize import fsolve    

class DAECalc:

    def __init__(self, probdef, dprobdef = None, d2probdef = None, doimport = True):

        if doimport:
            self.funcs = importlib.import_module(probdef)
        else:
            self.funcs = probdef
        self.var_types = list(self.funcs.var_types())
        self.res = self.funcs.res_function
        self.ic = self.funcs.ic_function
        self.jac = self.funcs.jac_function

        if dprobdef is not None:
            if doimport:
                self.dfuncs = importlib.import_module(dprobdef)
            else:
                self.dfuncs = dprobdef
            self.dres = self.dfuncs.res_function
            self.dic = self.dfuncs.ic_function
            self.djac = self.dfuncs.jac_function

        if d2probdef is not None:
            if doimport:
                self.d2funcs = importlib.import_module(d2probdef)
            else:
                self.d2funcs = d2probdef
            self.d2res = self.d2funcs.res_function
            self.d2ic = self.d2funcs.ic_function
            self.d2jac = self.d2funcs.jac_function
            
        
        self.atol = [1e-10]*len(self.var_types)
        self.rtol = [1e-10]*len(self.var_types)
        self.calculate_ic = True
        self.kwargs = {}

    def Calculateyprime(self, y, y_prime, pvalues, var_types, res, jac):

        def fitfunc(x, y, pvalues):
            y_prime = scipy.zeros(len(var_types))
            for i in range(len(var_types)):
                if var_types[i] == 1:
                    y_prime[i] = x[i]
                else:
                    y[i] = x[i]
            return res(0, y, y_prime, pvalues)


        def fitfuncjac(x, y, pvalues):
            y_prime = scipy.zeros(len(var_types))
            for i in range(len(var_types)):
                if var_types[i] == 1:
                    y_prime[i] = x[i]
                else:
                    y[i] = x[i]
            dy = jac(0, y, y_prime, 0.0, pvalues)
            dyprime = jac(0, y, y_prime, 1.0, pvalues) - dy
            ans = scipy.zeros((len(var_types),len(var_types)))
            for i in range(len(var_types)):
                if var_types[i] == 1:
                    ans[i] = dyprime[i]
                else:
                    ans[i] = dy[i]
            return ans


        x = y_prime*0 ## Initial guess

        for i in range(len(var_types)):
            if var_types[i] == 1:
                x[i] = y_prime[i]
            else:
                x[i] = y[i]

        x, info, ier, mesg = fsolve(fitfunc, x, (y, pvalues), fitfuncjac, full_output = 1, xtol = 1e-12) ## ideally we would like to specify an ftol
        
        for i in range(len(var_types)):
            if var_types[i] == 1:
                y_prime[i] = x[i]
            else:
                y[i] = x[i]

        return y, y_prime
                

    def evaluate(self,pvalues, ts, yp0=None):
        y0 = self.ic(pvalues)
        if yp0 == None:
            yp0 = 0.0*y0
        if self.calculate_ic == False:
            y0, yp0 = self.Calculateyprime(y0, yp0, pvalues, self.var_types, self.res, self.jac)
        return daeint(self.res, ts, y0, yp0, self.rtol, self.atol,
                      jac = self.jac, rpar = pvalues, var_types = self.var_types, calculate_ic = self.calculate_ic,
                      **self.kwargs)[0]

    def evaluate_derivative(self,pvalues, dpvalues, ts, yp0=None):
        y0 = self.dic(scipy.append(pvalues, dpvalues ))
        if yp0 == None:
            yp0 = y0*1.0
        if self.calculate_ic == False:
            y0, yp0 = self.Calculateyprime(y0, yp0, scipy.append(pvalues, dpvalues) , self.var_types*2, self.dres, self.djac)
        return daeint(self.dres, ts, y0, yp0, self.rtol*2, self.atol*2, calculate_ic = self.calculate_ic,
                      jac = self.djac, rpar = scipy.append(pvalues, dpvalues), var_types = self.var_types*2,
                      **self.kwargs)[0][:,len(y0)/2:]

    def evaluate_Avv(self,pvalues, dpvalues, ts, yp0=None):
        y0 = self.d2ic(list(pvalues) + list(dpvalues)*2 )
        if yp0 == None:
            yp0 = y0*0
        if self.calculate_ic == False:
            y0, yp0 = self.Calculateyprime(y0, yp0, list(pvalues) + list(dpvalues)*2, self.var_types*4, self.d2res, self.d2jac)
        return daeint(self.d2res, ts, y0, yp0, self.rtol*4, self.atol*4, calculate_ic = self.calculate_ic,
                      jac = self.d2jac, rpar = list(pvalues) + list(dpvalues)*2, var_types = self.var_types*4,
                      **self.kwargs)[0][:,3*len(y0)/4:]
    
    def evaluate_Auv(self,pvalues, dpvalues, dpvalues2, ts,yp0=None):
        y0 = self.d2ic(list(pvalues) + list(dpvalues) + list(dpvalues2) )
        if yp0 == None:
            yp0 = y0*0
        if self.calculate_ic == False:
            y0, yp0 = self.Calculateyprime(y0, yp0, list(pvalues) + list(dpvalues) + list(dpvalues2) , self.var_types*4, self.d2res, self.d2jac)
        return daeint(self.d2res, ts, y0, yp0, self.rtol*4, self.atol*4, calculate_ic = self.calculate_ic,
                      jac = self.d2jac, rpar = list(pvalues) + list(dpvalues) + list(dpvalues2), var_types = self.var_types*4,
                      **self.kwargs)[0][:,3*len(y0)/4:]
