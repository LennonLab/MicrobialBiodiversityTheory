from __future__ import division
from scipy import stats
import numpy as np

def e_simpson(SAD): # based on 1/D, not 1 - D

    " Simpson's evenness "
    SAD = filter(lambda a: a != 0, SAD)

    D = 0.0
    N = sum(SAD)
    S = len(SAD)

    for x in SAD:
        D += (x*x) / (N*N)

    E = round((1.0/D)/S, 4)

    if E < 0.0 or E > 1.0:
        print 'Simpsons Evenness =',E
    return E

def skewness(RAD):
    skew = stats.skew(RAD)
    # log-modulo skewness
    lms = np.log10(np.abs(float(skew)) + 1)
    if skew < 0:
        lms = lms * -1
    return lms
