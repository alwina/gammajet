import numpy as np
import scipy.integrate
import scipy.special

import iminuit


class FitFunction:
    def __init__(self):
        pass

    def getFunctionErrHigh(self):
        functionParams = {}

        for param in self.fitParams:
            functionParams[param] = self.fitParams[param] + self.fitErrors.get(param, 0)

        return self.getFunction(**functionParams)

    def getFunctionErrLow(self):
        functionParams = {}

        for param in self.fitParams:
            functionParams[param] = self.fitParams[param] - self.fitErrors.get(param, 0)

        return self.getFunction(**functionParams)

    def getFunctionFit(self):
        return self.getFunction(**self.fitParams)

    def isFiniteAndNonzero(self, y, yerr):
        finite = np.logical_and(np.isfinite(y), np.isfinite(yerr))
        return np.logical_and(finite, np.array(yerr) != 0)

    def calculateChi2(self, model, y, yerr):
        return np.sum(np.power(np.divide(np.subtract(y, model), yerr, where=self.isFiniteAndNonzero(y, yerr)), 2.0))

    def doFit(self, Chi2, nPoints):
        params = []
        for param, value in self.initParams.items():
            if not param.startswith('error_'):
                params.append(param)

        mt = iminuit.Minuit(Chi2, errordef=1, print_level=0, **self.initParams)
        mt.migrad()
        if not mt.migrad_ok():
            print('Warning: double ratio fit did not converge')

        self.fitParams = {}
        self.fitErrors = {}
        for param in params:
            self.fitParams[param] = mt.values[param]
            self.fitErrors[param] = mt.errors[param]

        self.chi2 = Chi2(**mt.values)
        self.dof = nPoints - len(params)
        self.chi2dof = self.chi2 / self.dof

        
class PromptPtFit(FitFunction):
    def __init__(self):
        self.initParams = {}
        self.initParams['A'] = 300
        self.initParams['error_A'] = 30
        self.initParams['p0'] = 100
        self.initParams['error_p0'] = 10
        self.initParams['p1'] = -100
        self.initParams['error_p1'] = 10
        self.initParams['p2'] = 100
        self.initParams['error_p2'] = 10
        self.initParams['p3'] = -10
        self.initParams['error_p3'] = 1
        self.initParams['p4'] = 1
        self.initParams['error_p4'] = 0.1
        
    def getFunction(self, A, p0, p1, p2, p3, p4):
        def function(x):
            return np.exp(p0 + p1 * np.log(x) + p2 * (np.log(x))**2 + p3 * (np.log(x))**3 + p4 * (np.log(x))**4) / A
        return function
    
    def getParamsAndErrors(self, y, yerr, x):
        def Chi2(A, p0, p1, p2, p3, p4):
            model = list(map(self.getFunction(A, p0, p1, p2, p3, p4), x))
            return self.calculateChi2(model, y, yerr)
        
        self.doFit(Chi2, len(y))


class NonPromptPtFit(FitFunction):
    def __init__(self):
        self.initParams = {}
        self.initParams['p0'] = -1e-05
        self.initParams['error_p0'] = 1e-06
        self.initParams['p1'] = 1e-03
        self.initParams['error_p1'] = 1e-04
        self.initParams['p2'] = -1e-02
        self.initParams['error_p2'] = 1e-03
        self.initParams['p3'] = 10
        self.initParams['error_p3'] = 1
        self.initParams['p4'] = 1
        self.initParams['error_p4'] = 0.1
        
    def getFunction(self, p0, p1, p2, p3, p4):
        prompt0 = 1.58972e+02
        prompt1 = -2.21542e+02
        prompt2 = 1.18650e+02
        prompt3 = -2.72830e+01
        prompt4 = 2.24137e+00
        promptA = 366.358683429
        
        def function(x):
            prompt = np.exp(prompt0 + prompt1 * np.log(x) + prompt2 * (np.log(x))**2 + prompt3 * (np.log(x))**3 + prompt4 * (np.log(x))**4) / promptA
            return prompt / (0.5 * (1 + scipy.special.erf(1e+4 * (p3 - x))) * (p0 * x**4 + p1 * x**3 + p2 * x**2 - 4 * (p0 * p3**3 + 3 * p1 * p3**2 + 2 * p2 * p3) * x + p4 + 3 * p0 * p3**4 + 2 * p1 * p3**3 + p2 * p3**2) + 0.5 * (1 + scipy.special.erf(1e+4 * (x - p3))) * p4)
        return function
    
    def getParamsAndErrors(self, y, yerr, x):
        def Chi2(p0, p1, p2, p3, p4):
            model = list(map(self.getFunction(p0, p1, p2, p3, p4), x))
            return self.calculateChi2(model, y, yerr)
        
        self.doFit(Chi2, len(y))
