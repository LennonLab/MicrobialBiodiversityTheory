from __future__ import division
import sys
import os
import numpy as np
from scipy import stats
import  matplotlib.pyplot as plt

import statsmodels.api as sm
import statsmodels.formula.api as smf
import pandas as pd
import time

########### PATHS ##############################################################

#tools = os.path.expanduser("~/tools")
# I (ken) have a tools folder in my home directory, where Cython compiled Python
# scripts for general 'tools' are stored. Keeps me from having a more than one
# copy of, say, macroeco_distributions lying around.
#sys.path.append(tools + "/macroeco_distributions")
#import macroeco_distributions as md

#macroeco = os.path.expanduser("~/GitHub/macroecotools")
#sys.path.append(macroeco + "/macroeco_distributions")
import macroeco_distributions as md
#sys.path.append(macroeco + "/macroecotools")
import macroecotools

######### CLASS ################################################################

class zipf:

    """ A class to obtain a zipf object with inherited mle shape parameter,
    mle form of the rank-abundance distribution, and a rank-abundance curve
    based on fitting the zipf to the observed via a generalized linear model."""

    def __init__(self, obs):
        self.obs = obs


    def from_cdf(self):
        """ Obtain the maximum likelihood form of the Zipf distribution, given
        the mle value for the Zipf shape parameter (a). Using a, this code
        generates a rank-abundance distribution (RAD) from the cumulative
        density function (cdf) using the percent point function (ppf) also known
        as the quantile function.
        see: http://www.esapubs.org/archive/ecol/E093/155/appendix-B.htm

        This is an actual form of the Zipf distribution, obtained from getting
        the mle for the shape parameter.
        """

        p = md.zipf_solver(self.obs)
        S = len(self.obs)
        rv = stats.zipf(a=p)
        rad = []

        for i in range(1, S+1):
            print rad
            val = (S - i + 0.5)/S
            x = rv.ppf(val)
            rad.append(int(x))


        return rad



    def from_glm(self):

        """ Fit the Zipf distribution to the observed vector of integer values
        using a generalized linear model.

        Note: This is a fitted curve; not an actual form of the Zipf distribution

        This method was inspired by the vegan
        package's open source code on vegan's public GitHub repository:
        https://github.com/vegandevs/vegan/blob/master/R/rad.zipf.R
        on Thursday, 19 Marth 2015 """

        ranks = np.log(range(1, len(self.obs)+1))
        off = [np.log(sum(self.obs))] * len(self.obs)

        d = pd.DataFrame({'ranks': ranks, 'off': off, 'x':self.obs})

        lm = smf.glm(formula='x ~ ranks', data = d, family = sm.families.Poisson()).fit()
        pred = lm.predict()

        return pred


ad = [20000, 10000, 8000, 6000, 1000, 200, 200, 100, 18, 16, 14, 12, 10, 4, 4, 2, 2, 2, 2, 2, 1,
            1, 1, 1, 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]

a = md.zipf_solver(ad)
S = len(ad)
rv = stats.zipf(a)


rad = []
vals = []
for i in range(1, S+1):
    vals.append((S - i + 0.5)/S)

t = time.clock()
x = rv.ppf(vals)
elapsed_t = time.clock() - t
print x, elapsed_t
#sys.exit()

ranks = range(1,len(ad)+1)

zipf_pred = zipf(ad)
zipf_cdf = zipf_pred.from_cdf()
#zipf_glm = zipf_pred.from_glm()
print zipf_pred
plt.scatter(ranks, np.log10(ad), color='blue') # observed rad
plt.scatter(ranks, np.log10(zipf_cdf), color='red') # mle cdf
#plt.scatter(ranks, np.log10(zipf_glm), color='green') # glm
#plt.plot(np.log10(zipf_cdf), np.log10(ad))
plt.savefig("zipf_test.png")

plt.show()
plt.close()
