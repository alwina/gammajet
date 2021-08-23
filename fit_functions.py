from __future__ import print_function

import numpy as np
import scipy.integrate
import scipy.misc
import scipy.special

import iminuit


def isFiniteAndNonzero(y, yerr):
    finite = np.logical_and(np.isfinite(y), np.isfinite(yerr))
    return np.logical_and(finite, np.array(yerr) != 0)


class FitFunction:
    def __init__(self):
        pass

    def getFunctionErrHigh(self):
        """
        Returns the function when varying all fit parameters up by their fit errors
        does not take into account any covariance
        """
        functionParams = {}

        for param in self.fitParams:
            functionParams[param] = self.fitParams[param] + self.fitErrors.get(param, 0)

        return self.getFunction(**functionParams)

    def getFunctionErrLow(self):
        """
        Returns the function when varying all fit parameters down by their fit errors
        does not take into account any covariance
        """
        functionParams = {}

        for param in self.fitParams:
            functionParams[param] = self.fitParams[param] - self.fitErrors.get(param, 0)

        return self.getFunction(**functionParams)

    def getFunctionFit(self):
        """Returns the function with the fit parameters"""
        return self.getFunction(**self.fitParams)

    def calculateChi2(self, model, y, yerr, xerr=None):
        """Generic function for calculating chi2"""
        if xerr is None:
            xerr = np.zeros_like(yerr)
        numerator = np.square(np.subtract(y, model))
        denominator = np.sum(np.square([yerr, xerr]), axis=0)
        divide = np.zeros_like(numerator)
        np.divide(numerator, denominator, out=divide, where=isFiniteAndNonzero(numerator, denominator))
        return np.sum(divide)

    def doFit(self, Chi2, nPoints):
        params = []
        for param, value in self.initParams.items():
            if not param.startswith('error_'):
                params.append(param)

        mt = iminuit.Minuit(Chi2, errordef=1, print_level=self.verbosity, **self.initParams)
        mt.migrad()
        if not mt.migrad_ok():
            print('Warning: fit did not converge')

        self.fitParams = {}
        self.fitErrors = {}
        for param in params:
            self.fitParams[param] = mt.values[param]
            self.fitErrors[param] = mt.errors[param]

        self.chi2 = Chi2(**mt.values)
        self.dof = nPoints - len(params)
        self.chi2dof = self.chi2 / self.dof


class SingleParameterLinearFit(FitFunction):
    def __init__(self, verbosity=0):
        self.initParams = {}
        self.initParams['m'] = 0.000001
        self.initParams['error_m'] = 0.1
        self.verbosity = verbosity

    def getFunction(self, m, x0, y0):
        def function(x):
            return m * (x - x0) + y0
        return function

    def getParamsAndErrors(self, y, yerr, x):
        x0 = np.mean([np.min(x), np.max(x)])
        weights = np.power(yerr, -2.0)
        weights[np.logical_not(np.isfinite(weights))] = 0
        y[np.logical_not(np.isfinite(y))] = 0
        y0 = np.average(y, weights=weights)

        def Chi2(m):
            model = list(map(self.getFunction(m, x0, y0), x))
            if self.verbosity == 1:
                print(self.calculateChi2(model, y, yerr))
            return self.calculateChi2(model, y, yerr)

        self.doFit(Chi2, len(y))

        self.fitParams['x0'] = x0
        self.fitParams['y0'] = y0

    def getResultText(self):
        return '$m={0:2.3f} \pm {1:2.3f}$'.format(self.fitParams['m'], self.fitErrors['m'])


class PowerLawFit(FitFunction):
    def __init__(self, verbosity=0):
        self.initParams = {}
        self.initParams['A'] = 10000
        self.initParams['error_A'] = 1000
        self.initParams['p'] = 4.0
        self.initParams['error_p'] = 0.4
        self.verbosity = verbosity

    def getFunction(self, A, p):
        def function(x):
            return A * np.power(x, -p)
        return function

    def getParamsAndErrors(self, y, yerr, x):
        def Chi2(A, p):
            model = list(map(self.getFunction(A, p), x))
            return self.calculateChi2(model, y, yerr)

        self.doFit(Chi2, len(y))


# must be normalized between xmin and xmax!
class SingleParameterPowerLawFit(FitFunction):
    def __init__(self, verbosity=0):
        self.initParams = {}
        self.initParams['p'] = 4.0
        self.initParams['error_p'] = 0.4
        self.verbosity = verbosity

    def getFunction(self, p, xmin, xmax):
        def function(x):
            integral = scipy.integrate.quad(lambda x: np.power(x, -p), xmin, xmax)[0]
            return np.power(x, -p) / integral
        return function

    def getParamsAndErrors(self, y, yerr, x, xmin, xmax):
        def Chi2(p):
            model = list(map(self.getFunction(p, xmin, xmax), x))
            return self.calculateChi2(model, y, yerr)

        self.doFit(Chi2, len(y))

        self.fitParams['xmin'] = xmin
        self.fitParams['xmax'] = xmax


class ErrorFunctionFit(FitFunction):
    def __init__(self, verbosity=0):
        self.initParams = {}
        self.initParams['a'] = 0.5
        self.initParams['error_a'] = 0.05
        self.initParams['b'] = 10.0
        self.initParams['error_b'] = 1.0
        self.initParams['c'] = 10.0
        self.initParams['error_c'] = 1.0
        self.verbosity = verbosity

    def getFunction(self, a, b, c):
        def function(x):
            return a * scipy.special.erf((x - b) / c)
        return function

    def getParamsAndErrors(self, y, yerr, x, xwidths=None):
        def Chi2(a, b, c):
            model = list(map(self.getFunction(a, b, c), x))
            if xwidths is None:
                xerr = np.zeros_like(x)
            else:
                xerr = []
                for (xval, xwidth) in zip(x, xwidths):
                    xerr.append(xwidth / 2 * scipy.misc.derivative(self.getFunction(a, b, c), xval))
            return self.calculateChi2(model, y, yerr, xerr)

        self.doFit(Chi2, len(y))
