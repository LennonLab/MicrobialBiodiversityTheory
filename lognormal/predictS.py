from __future__ import division
import numpy as np
import random
import scipy as sc
from scipy import stats
import os
import sys
import statsmodels.stats.api as sms
import statsmodels.formula.api as smf
from statsmodels.sandbox.regression.predstd import wls_prediction_std
from statsmodels.stats.outliers_influence import summary_table
from numpy import log, log2, exp, sqrt, log10, pi
from scipy.optimize import fsolve
import scipy.optimize as opt
import pandas as pd #import patsy
import mpmath as mpm
from scipy.optimize import fsolve
from math import erf, pi
import math

mydir = os.path.expanduser("~/GitHub/MicroMETE/")
mydir2 = os.path.expanduser("~/")
pi = math.pi

class predictS:
    def __init__(self, N, Nmax, predictNmax = True):
        self.N = N
        self.Nmax = Nmax
        self.predictNmax = predictNmax

    def alpha(self, a, Nmax, Nmin=1):

        """Numerically solve for Preston's a. Needed to estimate S using the lognormal"""

        y = sqrt(pi*Nmin*Nmax)/(2.0*a) * exp((a * log2(sqrt(Nmax/Nmin)))**2.0)
        y = y * exp((log(2.0)/(2.0*a))**2.0)
        y = y * erf(a * log2(sqrt(Nmax/Nmin)) - log(2.0)/(2.0*a))
        y += erf(a * log2(sqrt(Nmax/Nmin)) + log(2.0)/(2.0*a))
        y -= self.N

        return y # find alpha

    def s(self, a, Nmax, Nmin=1):

        """Predict S from the lognormal using Nmax as the only empirical input"""

        return sqrt(pi)/a * exp( (a * log2(sqrt(Nmax/Nmin)))**2) # Using equation 10


    def getNmax(self, b=0.6148, slope=0.942904468437):

        """Predict Nmax using N and the scaling law of Nmax with N predicted by the lognormal"""

        #NmaxCalc = 10 ** (b + slope*(log10(self.N)))
        NmaxCalc = (self.N **  slope) * b
        return int(round(NmaxCalc))


    def getS(self, predictNmax=True):

        guess = 0.1 # initial guess for Pre ston's alpha
        Nmin = 1

        if self.predictNmax == True:
            Nmax = self.getNmax()
        else:
            Nmax = self.Nmax

        a = opt.fsolve(self.alpha, guess, (Nmax, Nmin))[0]
        S2 = self.s(a, Nmax, 1)

        return int(round(S2))


#S = predictS(1000, 900, predictNmax=True).getS()
#print S
