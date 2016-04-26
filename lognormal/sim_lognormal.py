from __future__ import division
import numpy as np
import random
from random import randrange
import sys
import os

mydir = os.path.expanduser("~/")
sys.path.append(mydir + "GitHub/DiversityTools/macroeco_distributions")
sys.path.append(mydir + "GitHub/DiversityTools/macroecotools")
import macroeco_distributions
import macroecotools

from macroeco_distributions import pln, pln_solver


def get_rad_pln(S, mu, sigma, lower_trunc = True):
    """Obtain the predicted RAD from a Poisson lognormal distribution"""
    abundance = list(np.empty([S]))
    rank = range(1, int(S) + 1)
    cdf_obs = [(rank[i]-0.5) / S for i in range(0, int(S))]
    j = 0
    cdf_cum = 0
    i = 1
    while j < S:
        cdf_cum += pln.pmf(i, mu, sigma, lower_trunc)
        while cdf_cum >= cdf_obs[j]:
            abundance[j] = i
            j += 1
            if j == S:
                abundance.reverse()
                return abundance
        i += 1


def get_rad_from_obs(ab):
    mu, sigma = pln_solver(ab)
    pred_rad = get_rad_pln(len(ab), mu, sigma)
    return pred_rad




def AvgShape(RACs):
    """ Find the SAD in a random sample with the greatest average commonness
        among its ranked abundance states. This SAD is taken to represent the
        central tendency of the set, based on the SAD shape. """

    if len(RACs) > 1000:
        RACs = random.sample(RACs, 1000)

    a1 = 0 # SAD mean
    v1 = 0 # SAD variance
    for rac in RACs:
        in_common = []
        ct1 = 0
        for a in rac: # for each rank
            c = 0
            for sad in RACs:
                if a == sad[ct1]:
                    c += 1
            in_common.append(np.log(c))
            ct1 += 1
        a2 = np.mean(in_common)
        v2 = np.var(in_common)
        if a2 > a1:
            a1 = a2
            v1 = v2
            xRAD = rac
        elif a2 == a1:
            if v2 < v1:
                a1 = a2
                v1 = v2
                xRAD = rac

    return xRAD



def SimLogNormInt(N, S, sample_size):
    '''This script codes the Lognormal Model'''
    sample = []
    fails = 0

    while len(sample) < sample_size:

        n = int(round(0.75 * N))
        RAC = [n, N - n]

        while len(RAC) < S:

            ind = randrange(len(RAC))
            v = RAC.pop(ind)
            v1 = int(round(0.75 * v))
            v2 = v - v1   # forcing all abundance values to be integers

            if v1 < 1 or v2 < 1:
                fails += 1
                break  # forcing smallest abundance to be
                                        # greater than one
            RAC.extend([v1, v2])

        if len(RAC) == S and sum(RAC) == N:
            RAC.sort()
            RAC.reverse()
            sample.append(RAC)

    ad = AvgShape(sample)
    print len(sample),'SADs'

    return ad


def SimLogNormFloat(N, S, sample_size):
    '''This script codes the Lognormal Model'''
    sample = []
    while len(sample) < sample_size:

        n = 0.75 * N
        RAC = [n, N - n]

        while len(RAC) < S:

            ind = randrange(len(RAC))
            v = RAC.pop(ind)
            v1 = 0.75 * v
            v2 = v - v1   # forcing all abundance values to be integers

            RAC.extend([v1, v2])

        if len(RAC) == S and int(round(sum(RAC))) == N:
            RAC.sort()
            RAC.reverse()
            sample.append(RAC)

    ad = AvgShape(sample)

    return ad


"""
N = 1000
S = 50
sample_size = 1000

#SAD = SimLogNormFloat(N, S, sample_size)
SAD = SimLogNormInt(N, S, sample_size)
print 'N:', sum(SAD), 'S:', len(SAD), SAD

#plnSAD = get_rad_from_obs(SAD)
#print 'N:', sum(plnSAD), 'S:', len(plnSAD), plnSAD
"""
