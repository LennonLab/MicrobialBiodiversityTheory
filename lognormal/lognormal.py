from __future__ import division
import sys
import os
import numpy as np
from scipy import stats
import  matplotlib.pyplot as plt

from scipy.stats import norm

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

######### CLASS ################################################################

def ppoints(n):
    """ numpy analogue or `R`'s `ppoints` function
        see details at http://stat.ethz.ch/R-manual/R-patched/library/stats/html/ppoints.html
        :param n: array type or number

        Obtained from: http://stackoverflow.com/questions/20292216/imitating-ppoints-r-function-in-python
        on 5 April 2016
        """
    if n < 10:
        a = 3/8
    else:
        a = 1/2

    try:
        n = np.float(len(n))
    except TypeError:
        n = np.float(n)
    return (np.arange(n) + 1 - a)/(n + 1 - 2*a)



def lognorm_glm(ad):

    """ Fit the lognormal distribution to the observed vector of integer values
    using a generalized linear model.

    Note: This is a fitted curve; not an actual form of the lognormal distribution

    This method was inspired by the vegan package's open source code on vegan's public
    GitHub repository: https://github.com/vegandevs/vegan/R/rad.lognormal.R
    on Thursday, 5 April 2016
    """

    ranks = np.log(range(1, len(ad)+1))
    ranks = -norm.ppf(ppoints(len(ranks)))

    d = pd.DataFrame({'rnks': ranks, 'x': ad})
    lm = smf.glm(formula='x ~ rnks', data = d, family = sm.genmod.families.family.Poisson(link=sm.genmod.families.links.log)).fit()
    pred = lm.predict()

    return pred


ad = [20000, 10000, 8000, 6000, 1000, 200, 200, 100, 18, 16, 14, 12, 10, 4, 4, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]


ranks = range(1, len(ad)+1)
lnorm = lognorm_glm(ad)

plt.scatter(ranks, np.log10(ad), color='blue') # observed rad
plt.scatter(ranks, np.log10(lnorm), color='red') # mle cdf
#plt.savefig("lnorm_test.png")

plt.show()
